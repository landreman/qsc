# I don't really know what minimum version is required, but if you don't specify one, cmake gives a warning.
cmake_minimum_required (VERSION 3.9)

project(qsc LANGUAGES CXX)

find_package(MPI REQUIRED)

include_directories(include)

# file(GLOB SOURCES "src/*.cpp")

# Set where the compiled library will go.
# CMAKE_ARCHIVE_OUTPUT_DIRECTORY applies to static libraries.
# CMAKE_LIBRARY_OUTPUT_DIRECTORY applies to shared libraries.
# Set both:
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)

# Set where the executable will go:
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_library(qsc src/hw.cpp)
target_link_libraries(qsc PUBLIC MPI::MPI_CXX)

#add_executable(qsc_driver ${SOURCES})
add_executable(qsc_driver src/qsc.cpp)

# qsc driver depends on the qsc library:
target_link_libraries(qsc_driver PUBLIC qsc)
