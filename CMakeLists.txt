cmake_minimum_required(VERSION 3.22)



project(PyRasch C)

set(CMAKE_C_STANDARD 99)
include_directories(/opt/OpenBLAS/include)
link_directories(/opt/OpenBLAS/lib)

find_package(OpenMP REQUIRED)

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -march=native -mtune=native")

add_executable(PyRasch

        src/block.c
        src/block.h
        src/vector.c
        src/vector.h
        src/matrix.c
        src/matrix.h
        src/dichotomus.c
        src/dichotomus.h
        src/utils.c
        src/utils.h
        src/panopticon.c
        src/panopticon.h
        src/test.c)

target_link_libraries(PyRasch openblas m pthread OpenMP::OpenMP_C)

target_compile_options(PyRasch PRIVATE
        -march=native
        -mtune=native
        -flto
        -lpthread
        -fstrict-aliasing
        -funroll-loops
        -fopenmp
)