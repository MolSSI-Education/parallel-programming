cmake_minimum_required(VERSION 2.7)
project(MPI_EXAMPLE1)
enable_language(Fortran)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

add_executable(hello hello.cpp)
target_link_libraries(hello ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(hello PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(hello PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
