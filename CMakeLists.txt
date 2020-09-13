# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(gbwt)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

# Use openmp for parallelism, but it's configured differently on OSX
find_package(OpenMP)
if (OPENMP_FOUND)
  if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # assumes clang build
    # we can't reliably detect when we're using clang, so for the time being we assume
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
  endif()
endif()

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sdsl-lite"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

set(CMAKE_BUILD_TYPE Release)

add_library(gbwt STATIC
  ${CMAKE_SOURCE_DIR}/
  ${CMAKE_SOURCE_DIR}/algorithms.cpp
  ${CMAKE_SOURCE_DIR}/bwtmerge.cpp
  ${CMAKE_SOURCE_DIR}/cached_gbwt.cpp
  ${CMAKE_SOURCE_DIR}/dynamic_gbwt.cpp
  ${CMAKE_SOURCE_DIR}/fast_locate.cpp
  ${CMAKE_SOURCE_DIR}/files.cpp
  ${CMAKE_SOURCE_DIR}/gbwt.cpp
  ${CMAKE_SOURCE_DIR}/internal.cpp
  ${CMAKE_SOURCE_DIR}/metadata.cpp
  ${CMAKE_SOURCE_DIR}/support.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp
  ${CMAKE_SOURCE_DIR}/test.cpp
  ${CMAKE_SOURCE_DIR}/variants.cpp)
add_dependencies(gbwt sdsl-lite)

add_executable(build_gbwt ${CMAKE_SOURCE_DIR}/build_gbwt.cpp)
add_dependencies(build_gbwt gbwt)
target_include_directories(build_gbwt PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}")
target_link_libraries(build_gbwt
  "${LIBRARY_OUTPUT_PATH}/libgbwt.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a")

add_executable(build_ri ${CMAKE_SOURCE_DIR}/build_ri.cpp)
add_dependencies(build_ri gbwt)
target_include_directories(build_ri PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}")
target_link_libraries(build_ri
  "${LIBRARY_OUTPUT_PATH}/libgbwt.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a")

add_executable(merge_gbwt ${CMAKE_SOURCE_DIR}/merge_gbwt.cpp)
add_dependencies(merge_gbwt gbwt)
target_include_directories(merge_gbwt PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}")
target_link_libraries(merge_gbwt
  "${LIBRARY_OUTPUT_PATH}/libgbwt.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a")

add_executable(benchmark ${CMAKE_SOURCE_DIR}/benchmark.cpp)
add_dependencies(benchmark gbwt)
target_include_directories(benchmark PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}")
target_link_libraries(benchmark
  "${LIBRARY_OUTPUT_PATH}/libgbwt.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a")

add_executable(metadata_tool ${CMAKE_SOURCE_DIR}/metadata_tool.cpp)
add_dependencies(metadata_tool gbwt)
target_include_directories(metadata_tool PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}")
target_link_libraries(metadata_tool
  "${LIBRARY_OUTPUT_PATH}/libgbwt.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a")

add_executable(remove_seq ${CMAKE_SOURCE_DIR}/remove_seq.cpp)
add_dependencies(remove_seq gbwt)
target_include_directories(remove_seq PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}")
target_link_libraries(remove_seq
  "${LIBRARY_OUTPUT_PATH}/libgbwt.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a")

target_include_directories(gbwt PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}")
target_link_libraries(gbwt
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a")

if (APPLE)
elseif (TRUE)
  if (BUILD_STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
  endif()
endif()
