cmake_minimum_required(VERSION 3.22)

option(ASAN "Enable AddressSanitizer" OFF)
option(UBSAN "Enable UndefinedBehaviorSanitizer" OFF)


set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(elem_sym_poly C)

set(CMAKE_C_STANDARD 99)
include_directories(/opt/OpenBLAS/include ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/panopticon/include)
link_directories(/opt/OpenBLAS/lib ${CMAKE_CURRENT_SOURCE_DIR}/build/panopticon)

message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")

find_package(OpenMP REQUIRED)

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -march=native -mtune=native")

#add_subdirectory(panopticon)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)



add_library(elem_sym_poly SHARED
        src/polynomial.c
        #src/polynomial_test.c
        src/elementary_symmetric_polynomial.c
        src/mpfr_fft.c
        src/simd_fft.c
        src/elem_sym_poly.c
        include/elem_sym_poly.h
)

target_link_libraries(elem_sym_poly openblas m pthread OpenMP::OpenMP_C panopticon mpfr gmp)

set_target_properties(elem_sym_poly PROPERTIES PREFIX "")

target_compile_options(elem_sym_poly PRIVATE
        -O3
        -march=native
        -mtune=native
        -flto
        -lpthread
        -fopenmp
        -g
)

target_link_options(elem_sym_poly PRIVATE -flto)


if (ASAN)
    message("ADDRESS SANITIZER")
    target_compile_options(elem_sym_poly PRIVATE -g -fsanitize=address -fno-omit-frame-pointer)
    target_link_options(elem_sym_poly PRIVATE -fsanitize=address)
endif()

if (UBSAN)
    target_compile_options(elem_sym_poly PRIVATE  -g -fsanitize=undefined -fno-omit-frame-pointer)
    target_link_options(elem_sym_poly PRIVATE -fsanitize=undefined)
endif()