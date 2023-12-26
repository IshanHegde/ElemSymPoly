cmake_minimum_required(VERSION 3.4...3.18)

# All options are for development/testing purpose
option(ASAN "Enable AddressSanitizer" OFF)
option(UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(DEBUG "Debug mode" OFF)

project(pyElemSymPoly C)


set(CMAKE_C_STANDARD 99)

# For testing
if (DEBUG)
    SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
endif ()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

if (WIN32)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)
    pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
else()
    find_package(GMP REQUIRED)
    find_package(MPFR REQUIRED)
endif ()


include(FindPython)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${PYTHON_INCLUDE_DIR})

if(APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif()

add_library(pyElemSymPoly SHARED
        src/polynomial.c
        src/elementary_symmetric_polynomial.c
        src/mpfr_fft.c
        src/elem_sym_poly.c
)

if (APPLE)
    target_link_libraries(pyElemSymPoly m mpfr gmp ${Python_Library})
elseif(WIN32)
    target_link_libraries(pyElemSymPoly PkgConfig::mpfr PkgConfig::gmp ${Python_Library})
else()
    target_link_libraries(pyElemSymPoly m mpfr gmp ${Python_Library})
endif ()

set_target_properties(pyElemSymPoly PROPERTIES PREFIX "lib" )

target_compile_options(pyElemSymPoly PUBLIC
        -O3
        -flto
        -fomit-frame-pointer
        -fPIC
)

# Apple Platforms, by default, use the .dylib extension for shared libraries.
# CPython can only load .so files as modules on Linux and MacOS platforms.
if (APPLE)
    set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
    target_link_options(pyElemSymPoly PUBLIC -flto -fomit-frame-pointer -fPIC)
else()
    target_link_options(pyElemSymPoly PUBLIC -flto -fomit-frame-pointer -fPIC)
endif ()

target_compile_definitions(pyElemSymPoly PUBLIC VERSION_INFO=${EXAMPLE_VERSION_INFO})


if (ASAN)
    if(UNIX AND NOT APPLE)
        message("ADDRESS SANITIZER")
        target_compile_options(pyElemSymPoly PRIVATE -g -fsanitize=address -fno-omit-frame-pointer)
        target_link_options(pyElemSymPoly PRIVATE -fsanitize=address)
    else ()
        message("ADDRESS SANITIZER is currently only supported for Linux in this project")
    endif ()
endif()

if (UBSAN)
    if(UNIX AND NOT APPLE)
        message("UB SANITIZER")
        target_compile_options(pyElemSymPoly PRIVATE  -g -fsanitize=undefined -fno-omit-frame-pointer)
        target_link_options(pyElemSymPoly PRIVATE -fsanitize=undefined)
    else ()
        message("UB SANITIZER is currently only supported for Linux in this project")
    endif ()
endif()

#cibuildwheel --platform linux --archs all
#-DCMAKE_TOOLCHAIN_FILE=D:/a/ElemSymPoly/ElemSymPoly/vcpkg/scripts/buildsystems/vcpkg.cmake