cmake_minimum_required(VERSION 3.12)
project(omp_example3 CXX)

add_executable(example3 example3.cpp)

find_package(OpenMP REQUIRED)
target_link_libraries(example3 PUBLIC OpenMP::OpenMP_CXX)
