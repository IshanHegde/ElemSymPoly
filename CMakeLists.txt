cmake_minimum_required(VERSION 3.4...3.18)

option(ASAN "Enable AddressSanitizer" OFF)
option(UBSAN "Enable UndefinedBehaviorSanitizer" OFF)


#set(CMAKE_C_COMPILER "gcc")

project(pyElemSymPoly C)

set(CMAKE_C_STANDARD 99)

if (WIN32)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)
    pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
else()
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
    find_package(GMP REQUIRED)
    find_package(MPFR REQUIRED)
    find_package(Python3 COMPONENTS Development.Module REQUIRED)
endif ()



include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Python3_INCLUDE_DIRS})

message("Python_FOUND:${Python3_FOUND}")
message("Python_VERSION:${Python3_VERSION}")
message("Python_Development_FOUND:${Python3_Development_FOUND}")
message("Python_LIBRARIES:${Python3_LIBRARIES}")

if(APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
    include_directories(/usr/local/include)
    link_libraries(/usr/local/lib)
endif()

add_library(pyElemSymPoly SHARED
        src/polynomial.c
        src/elementary_symmetric_polynomial.c
        src/mpfr_fft.c
        src/elem_sym_poly.c
)

if (APPLE)
elseif(WIN32)
    target_link_libraries(pyElemSymPoly PkgConfig::mpfr PkgConfig::gmp ${Python3_LIBRARIES})
else()
    target_link_libraries(pyElemSymPoly m mpfr gmp)
endif ()

set_target_properties(pyElemSymPoly PROPERTIES PREFIX "lib" )

target_compile_options(pyElemSymPoly PUBLIC
        -O3
        -flto
        -fomit-frame-pointer
        -fPIC
)
if (APPLE)
    set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
else()
    target_link_options(pyElemSymPoly PUBLIC -flto -fomit-frame-pointer -fPIC)
endif ()

target_compile_definitions(pyElemSymPoly PUBLIC VERSION_INFO=${EXAMPLE_VERSION_INFO})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

if (ASAN)
    message("ADDRESS SANITIZER")
    target_compile_options(pyElemSymPoly PRIVATE -g -fsanitize=address -fno-omit-frame-pointer)
    target_link_options(pyElemSymPoly PRIVATE -fsanitize=address)
endif()

if (UBSAN)
    target_compile_options(pyElemSymPoly PRIVATE  -g -fsanitize=undefined -fno-omit-frame-pointer)
    target_link_options(pyElemSymPoly PRIVATE -fsanitize=undefined)
endif()

#cibuildwheel --platform linux --archs all
#-DCMAKE_TOOLCHAIN_FILE=D:/a/ElemSymPoly/ElemSymPoly/vcpkg/scripts/buildsystems/vcpkg.cmake