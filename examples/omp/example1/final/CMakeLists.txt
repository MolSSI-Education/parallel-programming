cmake_minimum_required(VERSION 3.12)
project(omp_example1 CXX)

add_executable(example1 example1.cpp)

find_package(OpenMP REQUIRED)
target_link_libraries(example1 PUBLIC OpenMP::OpenMP_CXX)
