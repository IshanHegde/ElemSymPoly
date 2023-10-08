cmake_minimum_required(VERSION 3.22)

option(ENABLE_ASAN "Enable AddressSanitizer" OFF)



set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(PyRasch C)

set(CMAKE_C_STANDARD 99)
include_directories(/opt/OpenBLAS/include ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/panopticon/include)
link_directories(/opt/OpenBLAS/lib ${CMAKE_CURRENT_SOURCE_DIR}/build/panopticon)

message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")

find_package(OpenMP REQUIRED)

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -march=native -mtune=native")

add_subdirectory(panopticon)

add_executable(PyRasch

        src/block.c
        src/vector.c
        src/matrix.c
        src/dichotomus.c
        src/utils.c
        src/temps.c
        #src/test_dft.c
        src/complex_block.c
        src/complex_vector.c
        src/roots.c
        src/dft.c
        src/arena.c)

target_link_libraries(PyRasch openblas m pthread OpenMP::OpenMP_C panopticon fftw3)

target_compile_options(PyRasch PRIVATE
        -O3
        -march=native
        -mtune=native
        -flto
        -lpthread
        -fstrict-aliasing
        -funroll-loops
        -fopenmp
        -ftree-vectorize
        -fomit-frame-pointer

)

if (ENABLE_ASAN)
    target_compile_options(PyRasch PRIVATE -fsanitize=address -fno-omit-frame-pointer)
    target_link_options(PyRasch PRIVATE -fsanitize=address)
endif()

