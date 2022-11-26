cmake_minimum_required (VERSION 3.12...3.18)
project(pyrasch)

include(ExternalProject)
set(EIGEN_BUILD_DIR   ${CMAKE_BINARY_DIR}/eigen)
set(EIGEN_INSTALL_DIR ${CMAKE_SOURCE_DIR}/include/eigen3)
set(EIGEN_INCLUDE_DIR ${EIGEN_INSTALL_DIR})

ExternalProject_Add(
    Eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
    PREFIX ${EIGEN_BUILD_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND
      ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/Eigen ${EIGEN_INCLUDE_DIR}/Eigen 
      && ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/unsupported ${EIGEN_INCLUDE_DIR}/unsupported
    TEST_COMMAND ""
)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})

ExternalProject_Add(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.2.1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)
include_directories(${source_dir}/include)



find_package(Eigen3 REQUIRED NO_MODULE)
find_package(pybind11 REQUIRED)

pybind11_add_module(pyrasch src/pyrasch.cpp)

add_dependencies(pyrasch pybind11)
add_dependencies(pyrasch Eigen)



target_compile_definitions(pyrasch PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  add_compile_options(-O3 -march=native -mtune=native)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-O3 -march=native -mtune=native)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  add_compile_options(-O3 -march=native -mtune=native)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  add_compile_options(-O3 -march=native -mtune=native)
endif()