cmake_minimum_required(VERSION 3.12)
project(omp_example2 CXX)

add_executable(example2 example2.cpp)

find_package(OpenMP REQUIRED)
target_link_libraries(example2 PUBLIC OpenMP::OpenMP_CXX)
