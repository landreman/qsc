# I don't really know what minimum version is required, but if you don't specify one, cmake gives a warning.
# Versions >= 3.9 have improved support for MPI: https://cliutils.gitlab.io/modern-cmake/chapters/packages/MPI.html
cmake_minimum_required (VERSION 3.9)

# Comment out the next line when speed is needed
set(CMAKE_BUILD_TYPE Debug)

# Figure out the command to submit mpi jobs. It should be either srun or mpiexec. Other batch systems besides Slurm are not supported.
# First, see if the system uses slurm by checking for the "sinfo" command.
execute_process(COMMAND "sinfo" RESULT_VARIABLE SINFO_RESULT OUTPUT_QUIET ERROR_QUIET)
message("Result of 'sinfo': ${SINFO_RESULT}")
# If sinfo ran successfully, the "result_variable" SINFO_RESULT will be 0.
# Next, see if mpiexec exists:
execute_process(COMMAND mpiexec --version RESULT_VARIABLE MPIEXEC_RESULT OUTPUT_QUIET ERROR_QUIET)
message("Result of 'mpiexec --version': ${MPIEXEC_RESULT}")
# If mpiexec ran successfully, the "result_variable" MPIEXEC_RESULT will be 0.
if (SINFO_RESULT EQUAL 0)
  message("srun detected")
  set(QSC_COMMAND_TO_SUBMIT_JOB "srun -n NUM_PROCS")
elseif (MPIEXEC_RESULT EQUAL 0)
  message("srun not detected. mpiexec detected")
  set(QSC_COMMAND_TO_SUBMIT_JOB "mpiexec -n NUM_PROCS")
else()
  message("WARNING!! Neither slurm nor mpiexec was detected, so I am not sure how to launch MPI jobs.")
endif()

#set(CMAKE_MODULE_PATH cmake)
include(cmake/knownHosts.cmake)

message("QSC_COMMAND_TO_SUBMIT_JOB is ${QSC_COMMAND_TO_SUBMIT_JOB}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/tests/commandToSubmitJob ${QSC_COMMAND_TO_SUBMIT_JOB})

project(qsc LANGUAGES CXX)

message("BLA_VENDOR: ${BLA_VENDOR}")
find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)

# Tell "make" to print out the commands used for compiling and linking:
set(CMAKE_VERBOSE_MAKEFILE on)

include_directories(include)
include_directories(src)
include_directories(externalPackages/doctest)

# Get a list of all source files except for the driver:
file(GLOB SOURCES "src/*.cpp")
message("Sources before filter: ${SOURCES}")
list(FILTER SOURCES EXCLUDE REGEX ".*qsc_driver.cpp$")
message("Sources: ${SOURCES}")

# Set where the compiled library will go.
# CMAKE_ARCHIVE_OUTPUT_DIRECTORY applies to static libraries.
# CMAKE_LIBRARY_OUTPUT_DIRECTORY applies to shared libraries.
# Set both:
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)

# Set where the executable will go:
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_library(qsc ${SOURCES})
# Below, PUBLIC means that anything that links to qsc must also link to MPI.
target_link_libraries(qsc PUBLIC MPI::MPI_CXX ${BLAS_LIBRARIES})

add_executable(qsc_driver src/qsc_driver.cpp)

# qsc driver depends on the qsc library:
target_link_libraries(qsc_driver PUBLIC qsc)

# Set up unitTests executable
file(GLOB UNIT_TEST_SOURCES "src/tests/*.cpp")
message("Tests: ${UNIT_TEST_SOURCES}")
add_executable(unitTests ${UNIT_TEST_SOURCES})
target_link_libraries(unitTests PUBLIC qsc)
# Put the unitTests executable in the "tests" directory:
set_property(TARGET unitTests PROPERTY RUNTIME_OUTPUT_DIRECTORY tests)
# Doctest requires c++11
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

# Copy some files from the source to the build locations
message("CMAKE_CURRENT_SOURCE_DIR is ${CMAKE_CURRENT_SOURCE_DIR}")
message("CMAKE_CURRENT_BINARY_DIR is ${CMAKE_CURRENT_BINARY_DIR}")
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message("Building in root directory, so not copying files from source to build location.")
else()
  message("Copying files from source to build location.")
  message("CMake version: ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")
  # It is preferable to create links, but CMake only introduced links in version 3.14
  if ((CMAKE_MAJOR_VERSION GREATER 3) OR (CMAKE_MINOR_VERSION GREATER 13))
    message("  using links")
    file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/tests/runTests ${CMAKE_CURRENT_BINARY_DIR}/tests/runTests COPY_ON_ERROR SYMBOLIC)
    file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/tests/runMPIUnitTests ${CMAKE_CURRENT_BINARY_DIR}/tests/runMPIUnitTests COPY_ON_ERROR SYMBOLIC)
  else()
    message("  using copy")
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests/runTests DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests/runMPIUnitTests DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests)
  endif()
endif()


find_package(GSL)
message("GSL_FOUND = ${GSL_FOUND}")
if(GSL_FOUND)
  message("Found GSL")
else()
  message("Could not find GSL")
endif()

