cmake_minimum_required(VERSION 3.13)
PROJECT(pathIntegral CXX)

# settings to add -std=c++17 option, and avoid -std=gnu++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# settings for fortran
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2018")
SET (CMAKE_Fortran_COMPILER  gfortran)
ENABLE_LANGUAGE(Fortran)

set(CMAKE_EXPORT_COMPILE_COMMANDS OFF) # ON to get a .json with compiler calls for all targets

find_library(GMSH_LIB gmsh)
if(GMSH_LIB)
    message(STATUS "Found GMSH")
else()
    message(FATAL_ERROR "Could not find libgmsh")
endif()

find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()

#===========================================================================
# Option to enable AddressSanitizer
# cmake . -DENABLE_SANITIZER=ON
#===========================================================================
option(ENABLE_SANITIZER "Enable AddressSanitizer" OFF)
message(STATUS "Enable AddressSanitizer: ENABLE_SANITIZER = ${ENABLE_SANITIZER}")

if(ENABLE_SANITIZER)
    set(SANITIZER_FLAGS "-fsanitize=address")
    add_compile_options(${SANITIZER_FLAGS})
    add_link_options(${SANITIZER_FLAGS})
endif()

#===========================================================================
# cmake . -DENABLE_SEQ=OFF
# execution policy, by default it is parallel, sequential is for debug purpose.
# when ENABLE_SEQ is ON, executable does not link to TBB
#===========================================================================

option(ENABLE_SEQ "Enable sequential execution policy" OFF)
message( STATUS "Enable sequential execution policy: ENABLE_SEQ = " ${ENABLE_SEQ} )

if(NOT ${ENABLE_SEQ})
    find_package(TBB QUIET)
    if(TBB_FOUND)
        message(STATUS "Found TBB, parallel mode is active")
        set(TBB_LINK TBB::tbb)
    else()
        message(STATUS "TBB NOT FOUND")
        message(WARNING "sequential mode is active")
        set(ENABLE_SEQ TRUE)
    endif()
endif()

#===========================================================================
# TinyXML2 integration (installed in /usr/local)
#===========================================================================
# Add /usr/local to the search paths
set(CMAKE_PREFIX_PATH "/usr/local" ${CMAKE_PREFIX_PATH})
set(CMAKE_LIBRARY_PATH "/usr/local/lib" ${CMAKE_LIBRARY_PATH})
set(CMAKE_INCLUDE_PATH "/usr/local/include" ${CMAKE_INCLUDE_PATH})

# Locate tinyxml2 manually
find_path(TINYXML2_INCLUDE_DIR tinyxml2.h PATHS /usr/local/include)
if(NOT TINYXML2_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find tinyxml2.h in /usr/local/include")
else()
    message(STATUS "Found tinyxml2 include directory: ${TINYXML2_INCLUDE_DIR}")
    include_directories(${TINYXML2_INCLUDE_DIR})
endif()

find_library(TINYXML2_LIBRARY tinyxml2 PATHS /usr/local/lib)
if(NOT TINYXML2_LIBRARY)
    message(FATAL_ERROR "Could not find the tinyxml2 library in /usr/local/lib")
else()
    message(STATUS "Found tinyxml2 library: ${TINYXML2_LIBRARY}")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# search for pkg-config
include (FindPkgConfig)
if (NOT PKG_CONFIG_FOUND)
    message (FATAL_ERROR "pkg-config not found")
endif ()

# check for libpng
pkg_check_modules (LIBPNG libpng16 REQUIRED)
if (NOT LIBPNG_FOUND)
    message(FATAL_ERROR "You don't seem to have libpng16 development libraries installed")
else ()
    include_directories (${LIBPNG_INCLUDE_DIRS})
    link_directories (${LIBPNG_LIBRARY_DIRS})
    link_libraries (${LIBPNG_LIBRARIES})
endif ()

#===========================================================================
# these parameters are used to build a config.h file 
#===========================================================================

execute_process(
    COMMAND git describe --tags --always --dirty
    OUTPUT_VARIABLE PATHINTEGRAL_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE
)
# remove leading "v" from tag name
string(REGEX REPLACE "^v" "" PATHINTEGRAL_VERSION ${PATHINTEGRAL_VERSION})

message(STATUS "pathIntegral version = ${PATHINTEGRAL_VERSION}")

# Link-time optimization requires cmake version at least 3.9.
option(ENABLE_LTO "Enable link-time optimization" OFF)
message(STATUS "Link-time optimization: ENABLE_LTO = " ${ENABLE_LTO})
if(${ENABLE_LTO})
    cmake_policy(SET CMP0069 NEW)
    include(CheckIPOSupported)
    check_ipo_supported()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

SET(HEADERS config.h node.h  expression_parser.h mesh.h
    termcolor.h settings.h tetra.h facette.h tags.h
    chronometer.h element.h geo.h fem2d.h write_png.h pathIntegral.h)

SET(SOURCES settings.cpp read.cpp tetra.cpp expression_parser.cpp chronometer.cpp
    tags.cpp write_png.cpp exportRatioGrayScaleImage.cpp exportMagIntegrals.cpp chapeaux.cpp fmm2d_sum.cpp femutil.cpp correction.cpp 
    exportHoloPhase.cpp exportRatioRGBscaleImage.cpp pathIntegral.cpp)

set(SOURCES_FORTRAN rfmm2d.f90)

configure_file(config.h.in ./config.h)

add_custom_command(
    OUTPUT default-settings.o
    COMMAND ld -r -b binary -z noexecstack default-settings.yml -o default-settings.o
    DEPENDS default-settings.yml)
    
#===========================================================================
# Compiler options
#===========================================================================
# -ffast-math -msse4.2 -mfpmath=sse (on X86 -mfpmath=sse set by default) or -mavx
# -fsanitize=address leak null undefined # both compiler and linker
#===========================================================================

ADD_EXECUTABLE(pathIntegral main.cpp ${HEADERS} ${SOURCES} ${SOURCES_FORTRAN} default-settings.o)
target_compile_options( pathIntegral PUBLIC -Wall -Wextra -O3 -march=native )
target_compile_definitions(pathIntegral PUBLIC -DCGAL_HAS_THREADS -DNDEBUG -DTBB_SUPPRESS_DEPRECATED_MESSAGES )
target_link_libraries(pathIntegral yaml-cpp duktape ${TBB_LINK} Eigen3::Eigen ${GMSH_LIB} ${TINYXML2_LIBRARY} 
	-lfmm2d -lgfortran -lm -lblas -llapack)

#===========================================================================

install(TARGETS pathIntegral RUNTIME DESTINATION bin)

