/usr/bin/cmake -S/home/toussaij/Devel/Diwina_path_integral -B/home/toussaij/Devel/Diwina_path_integral --check-build-system CMakeFiles/Makefile.cmake 0
/usr/bin/cmake -E cmake_progress_start /home/toussaij/Devel/Diwina_path_integral/CMakeFiles /home/toussaij/Devel/Diwina_path_integral//CMakeFiles/progress.marks
make  -f CMakeFiles/Makefile2 all
make[1]: Entering directory '/home/toussaij/Devel/Diwina_path_integral'
make  -f CMakeFiles/pathIntegral.dir/build.make CMakeFiles/pathIntegral.dir/depend
make[2]: Entering directory '/home/toussaij/Devel/Diwina_path_integral'
[  4%] Generating default-settings.o
ld -r -b binary -z noexecstack default-settings.yml -o default-settings.o
cd /home/toussaij/Devel/Diwina_path_integral && /usr/bin/cmake -E cmake_depends "Unix Makefiles" /home/toussaij/Devel/Diwina_path_integral /home/toussaij/Devel/Diwina_path_integral /home/toussaij/Devel/Diwina_path_integral /home/toussaij/Devel/Diwina_path_integral /home/toussaij/Devel/Diwina_path_integral/CMakeFiles/pathIntegral.dir/DependInfo.cmake "--color="
Dependee "/home/toussaij/Devel/Diwina_path_integral/CMakeFiles/pathIntegral.dir/DependInfo.cmake" is newer than depender "/home/toussaij/Devel/Diwina_path_integral/CMakeFiles/pathIntegral.dir/depend.internal".
Dependee "/home/toussaij/Devel/Diwina_path_integral/CMakeFiles/CMakeDirectoryInformation.cmake" is newer than depender "/home/toussaij/Devel/Diwina_path_integral/CMakeFiles/pathIntegral.dir/depend.internal".
Scanning dependencies of target pathIntegral
make[2]: Leaving directory '/home/toussaij/Devel/Diwina_path_integral'
make  -f CMakeFiles/pathIntegral.dir/build.make CMakeFiles/pathIntegral.dir/build
make[2]: Entering directory '/home/toussaij/Devel/Diwina_path_integral'
[  9%] Building Fortran object CMakeFiles/pathIntegral.dir/rfmm2d.f90.o
/usr/bin/gfortran -DCGAL_HAS_THREADS -DCPP17 -DENABLE_XTL_COMPLEX -DHAVE_CBLAS=1 -DNDEBUG -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DXSIMD_ENABLE_XTL_COMPLEX -DXTENSOR_DISABLE_EXCEPTIONS -DXTENSOR_FFTW_USE_DOUBLE -DXTENSOR_FFTW_USE_FLOAT -DXTENSOR_USE_XSIMD -I/home/toussaij/Devel/Diwina_path_integral/.. -I/toto -I/home/toussaij/Devel/ScalFMM/include -I/home/toussaij/Devel/ScalFMM/modules/internal/inria_tools -I/home/toussaij/Devel/ScalFMM/modules/internal/xtensor-fftw/include -I/home/toussaij/Devel/ScalFMM/modules/internal/cpp_tools/colors/include -I/home/toussaij/Devel/ScalFMM/modules/internal/cpp_tools/timers/include -I/home/toussaij/Devel/ScalFMM/modules/internal/cpp_tools/cl_parser/include -I/home/toussaij/Devel/ScalFMM/modules/internal/xsimd/include -I/home/toussaij/Devel/ScalFMM/modules/internal/xtensor/include -I/home/toussaij/Devel/ScalFMM/modules/internal/xtensor-blas/include -I/home/toussaij/Devel/ScalFMM/modules/internal/xtl/include -I/usr/include/libpng16 -I/usr/include -I/usr/include/eigen3 -std=f2018 -Wall -Wextra -O3 -march=native -c /home/toussaij/Devel/Diwina_path_integral/rfmm2d.f90 -o CMakeFiles/pathIntegral.dir/rfmm2d.f90.o
f951: Warning: Nonexistent include directory ‘/toto’ [-Wmissing-include-dirs]
/usr/bin/cmake -E cmake_copy_f90_mod glue.mod CMakeFiles/pathIntegral.dir/glue.mod.stamp GNU
/usr/bin/cmake -E touch CMakeFiles/pathIntegral.dir/rfmm2d.f90.o.provides.build
[ 13%] Building CXX object CMakeFiles/pathIntegral.dir/main.cpp.o
x86_64-w64-mingw32-g++ -DCGAL_HAS_THREADS -DCPP17 -DENABLE_XTL_COMPLEX -DHAVE_CBLAS=1 -DNDEBUG -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DXSIMD_ENABLE_XTL_COMPLEX -DXTENSOR_DISABLE_EXCEPTIONS -DXTENSOR_FFTW_USE_DOUBLE -DXTENSOR_FFTW_USE_FLOAT -DXTENSOR_USE_XSIMD -I/home/toussaij/Devel/Diwina_path_integral/.. -I/toto -I/home/toussaij/Devel/ScalFMM/include -I/home/toussaij/Devel/ScalFMM/modules/internal/inria_tools -I/home/toussaij/Devel/ScalFMM/modules/internal/xtensor-fftw/include -I/home/toussaij/Devel/ScalFMM/modules/internal/cpp_tools/colors/include -I/home/toussaij/Devel/ScalFMM/modules/internal/cpp_tools/timers/include -I/home/toussaij/Devel/ScalFMM/modules/internal/cpp_tools/cl_parser/include -I/home/toussaij/Devel/ScalFMM/modules/internal/xsimd/include -I/home/toussaij/Devel/ScalFMM/modules/internal/xtensor/include -I/home/toussaij/Devel/ScalFMM/modules/internal/xtensor-blas/include -I/home/toussaij/Devel/ScalFMM/modules/internal/xtl/include -I/usr/include/libpng16 -isystem /usr/include/eigen3 -fopenmp -static-libstdc++ -std=c++17 -Wall -Wextra -O3 -march=native -MD -MT CMakeFiles/pathIntegral.dir/main.cpp.o -MF CMakeFiles/pathIntegral.dir/main.cpp.o.d -o CMakeFiles/pathIntegral.dir/main.cpp.o -c /home/toussaij/Devel/Diwina_path_integral/main.cpp
In file included from /home/toussaij/Devel/Diwina_path_integral/mesh.h:17,
                 from /home/toussaij/Devel/Diwina_path_integral/fem.h:25,
                 from /home/toussaij/Devel/Diwina_path_integral/main.cpp:10:
/home/toussaij/Devel/Diwina_path_integral/geo.h:1:10: fatal error: CGAL/Simple_cartesian.h: No such file or directory
    1 | #include <CGAL/Simple_cartesian.h>
      |          ^~~~~~~~~~~~~~~~~~~~~~~~~
compilation terminated.
make[2]: *** [CMakeFiles/pathIntegral.dir/build.make:80: CMakeFiles/pathIntegral.dir/main.cpp.o] Error 1
make[2]: Leaving directory '/home/toussaij/Devel/Diwina_path_integral'
make[1]: *** [CMakeFiles/Makefile2:83: CMakeFiles/pathIntegral.dir/all] Error 2
make[1]: Leaving directory '/home/toussaij/Devel/Diwina_path_integral'
make: *** [Makefile:136: all] Error 2
