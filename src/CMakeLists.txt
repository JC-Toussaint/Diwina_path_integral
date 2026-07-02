SET(HEADERS config.h node.h expression_parser.h mesh.h
    termcolor.h settings.h tetra.h facette.h tags.h
    chronometer.h element.h geo.h fem2d.h write_png.h pathIntegral.h
    scalfmm_2d.hpp scalfmm3.hpp )

SET(SOURCES settings.cpp read.cpp tetra.cpp expression_parser.cpp chronometer.cpp
    tags.cpp write_png.cpp exportRatioGrayScaleImage.cpp exportMagIntegrals.cpp chapeaux.cpp
    scalfmm2d_sum.cpp direct2d_sum.cpp femutil.cpp correction.cpp exportHoloPhase.cpp
    exportMagIntegrals_nd_HoloPhase.cpp pathIntegral.cpp)

configure_file(config.h.in ./config.h)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

add_custom_command(
    OUTPUT default-settings.o
    COMMAND ld -r -b binary -z noexecstack default-settings.yml -o default-settings.o
    DEPENDS default-settings.yml)

include_directories(${CMAKE_SOURCE_DIR})

ADD_EXECUTABLE(pathIntegral main.cpp ${HEADERS} ${SOURCES} ${SOURCES_FORTRAN} default-settings.o)
target_compile_options( pathIntegral PUBLIC -Wall -Wextra -O3 -march=native )
target_compile_definitions(pathIntegral PUBLIC CGAL_HAS_THREADS NDEBUG TBB_SUPPRESS_DEPRECATED_MESSAGES)

# Link - order is important to avoid dependencies problems
target_link_libraries(pathIntegral
    yaml-cpp
    duktape
    ${TBB_LINK}
    Eigen3::Eigen
    gmsh_target
    scalfmm
    LAPACK::LAPACK
    BLAS::BLAS
    ${FFTW_TARGET}
    OpenMP::OpenMP_CXX
    -ltbb
    -lm
    -lboost_program_options
)

