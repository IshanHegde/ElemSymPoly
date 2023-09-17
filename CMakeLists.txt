cmake_minimum_required(VERSION 3.26)



project(item C)

set(CMAKE_C_STANDARD 99)
include_directories(/opt/OpenBLAS/include)
link_directories(/opt/OpenBLAS/lib)

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -march=native -mtune=native")

add_executable(item
        src/matrix.c
        src/matrix.h
        src/test.c)

target_link_libraries(item openblas)

target_compile_options(item PRIVATE
        -march=native
        -mtune=native
        -flto
        -lpthread
        -fstrict-aliasing
        -funroll-loops
)