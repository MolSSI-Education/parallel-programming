cmake_minimum_required(VERSION 2.7)
project(MPI_EXAMPLE3)
enable_language(Fortran)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

add_executable(mc mc.cpp)
target_link_libraries(mc ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(mc PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(mc PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
