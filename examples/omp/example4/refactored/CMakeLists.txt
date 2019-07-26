cmake_minimum_required(VERSION 3.12)
project(omp_example4 CXX)

add_executable(example4 example4.cpp)

find_package(OpenMP REQUIRED)
target_link_libraries(example4 PUBLIC OpenMP::OpenMP_CXX)
