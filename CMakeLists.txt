cmake_minimum_required(VERSION 3.6)
project(abr_omp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES abr-omp.cpp)
add_executable(abr_omp ${SOURCE_FILES})