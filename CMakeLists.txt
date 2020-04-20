# I don't really know what minimum version is required, but if you don't specify one, cmake gives a warning.
# Versions >= 3.9 have improved support for MPI: https://cliutils.gitlab.io/modern-cmake/chapters/packages/MPI.html
cmake_minimum_required (VERSION 3.9)

project(qsc LANGUAGES CXX)

find_package(MPI REQUIRED)

set(CMAKE_VERBOSE_MAKEFILE on)

include_directories(include)
include_directories(externalPackages/catch2)

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
# Below, PUBLIC means that anything that links to qsc must also link to MPI.
target_link_libraries(qsc PUBLIC MPI::MPI_CXX)

#add_executable(qsc_driver ${SOURCES})
add_executable(qsc_driver src/qsc.cpp)

# qsc driver depends on the qsc library:
target_link_libraries(qsc_driver PUBLIC qsc)

# Set up unitTests executable
file(GLOB UNIT_TEST_SOURCES "src/tests/*.cpp")
add_executable(unitTests ${UNIT_TEST_SOURCES})
target_link_libraries(unitTests PUBLIC qsc)
# Put the unitTests executable in the "tests" directory:
set_property(TARGET unitTests PROPERTY RUNTIME_OUTPUT_DIRECTORY tests)
# Catch2 requires c++11
set_property(TARGET unitTests PROPERTY CXX_STANDARD 11)
# Do not build the unit tests by default, since they take some time to build:
set_property(TARGET unitTests PROPERTY EXCLUDE_FROM_ALL 1)



#enable_testing()
#
##add_test(unitTests unitTests)
#add_test(NAME unitTests COMMAND unitTests)
##add_test(mpiUnitTests tests/runMPIUnitTests)
#add_test(NAME mpiUnitTests COMMAND runMPIUnitTests WORKING_DIRECTORY tests)

# Set up "make test" so it runs the scripts I want it to run.
# See https://cmake.org/cmake/help/latest/command/add_custom_target.html
add_custom_target(test COMMAND ./runTests WORKING_DIRECTORY tests DEPENDS unitTests)
