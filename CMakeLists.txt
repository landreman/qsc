# I don't really know what minimum version is required, but if you don't specify one, cmake gives a warning.
# Versions >= 3.9 have improved support for MPI: https://cliutils.gitlab.io/modern-cmake/chapters/packages/MPI.html
# The command add_compile_definitions, used in this file, was new in CMake version 3.12.
cmake_minimum_required (VERSION 3.12)

# Comment out the next line when speed is needed
#set(CMAKE_BUILD_TYPE Debug)
add_definitions(-O2)
#add_definitions(-fsanitize=address -fno-omit-frame-pointer)

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

# If on a Mac, use the Accelerate framework for BLAS/LAPACK.
# Otherwise, CMake may use a random BLAS/LAPACK it finds that is not as optimized.
# To check if we're on a Mac, evaluate 'uname -s'.
execute_process(COMMAND uname -s
                OUTPUT_VARIABLE UNAME_S
                OUTPUT_STRIP_TRAILING_WHITESPACE)
message("Result of uname -s: ${UNAME_S}")
if ("${UNAME_S}" STREQUAL "Darwin")
  message("Setting BLA_VENDOR to Apple")
  set(BLA_VENDOR Apple)
endif()

#set(CMAKE_MODULE_PATH cmake)
include(cmake/knownHosts.cmake)

project(qsc LANGUAGES CXX)

#set(CMAKE_MODULE_PATH cmake)
# Tell CMake where to find FindNetCDF.cmake:
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#set (NETCDF_CXX "YES")
find_package (NetCDF REQUIRED)
message("NETCDF_INCLUDES=${NETCDF_INCLUDES}")
message("NETCDF_LIBRARIES=${NETCDF_LIBRARIES}")

message("QSC_COMMAND_TO_SUBMIT_JOB is ${QSC_COMMAND_TO_SUBMIT_JOB}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/tests/commandToSubmitJob ${QSC_COMMAND_TO_SUBMIT_JOB})

message("BLA_VENDOR: ${BLA_VENDOR}")
find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GSL REQUIRED)

# Tell "make" to print out the commands used for compiling and linking:
set(CMAKE_VERBOSE_MAKEFILE on)

include_directories(include)
include_directories(src)
include_directories(externalPackages/doctest)
include_directories(externalPackages/toml11)
include_directories("${NETCDF_INCLUDES}")
include_directories("${GSL_INCLUDE_DIRS}")

# Get a list of all source files except for the driver:
file(GLOB SOURCES "src/*.cpp")
message("Sources before filter: ${SOURCES}")
list(FILTER SOURCES EXCLUDE REGEX ".*xqsc.cpp$")
message("Sources: ${SOURCES}")

# Set where the compiled library will go.
# CMAKE_ARCHIVE_OUTPUT_DIRECTORY applies to static libraries.
# CMAKE_LIBRARY_OUTPUT_DIRECTORY applies to shared libraries.
# Set both:
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)

# Set where the executable will go:
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# single precision is off by default
option(SINGLE "Compile in single precision" OFF)
if (SINGLE)
  message("Compiling in SINGLE precision")
  # Add -DFOO to the compile statements
  add_compile_definitions(SINGLE)
  set(QSC_LIB qsc_single)
  set(QSC_DRIVER xqsc_single)
  set(QSC_TESTS unitTests_single)
else()
  message("Compiling in DOUBLE precision")
  set(QSC_LIB qsc)
  set(QSC_DRIVER xqsc)
  set(QSC_TESTS unitTests)
endif()

add_library(${QSC_LIB} ${SOURCES})
# Below, PUBLIC means that anything that links to qsc must also link to MPI, BLAS, & LAPACK.
target_link_libraries(${QSC_LIB} PUBLIC MPI::MPI_CXX ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${GSL_LIBRARIES})

# toml11 requires c++11
set_property(TARGET ${QSC_LIB} PROPERTY CXX_STANDARD 11)

add_executable(${QSC_DRIVER} src/xqsc.cpp)
set_property(TARGET ${QSC_DRIVER} PROPERTY CXX_STANDARD 11)

# qsc driver depends on the qsc library:
target_link_libraries(${QSC_DRIVER} PUBLIC ${QSC_LIB})

# Set up unitTests executable
file(GLOB UNIT_TEST_SOURCES "src/tests/*.cpp")
message("Tests: ${UNIT_TEST_SOURCES}")
add_executable(${QSC_TESTS} ${UNIT_TEST_SOURCES})
target_link_libraries(${QSC_TESTS} PUBLIC ${QSC_LIB})
# Put the unitTests executable in the "tests" directory:
set_property(TARGET ${QSC_TESTS} PROPERTY RUNTIME_OUTPUT_DIRECTORY tests)
# Doctest requires c++11
set_property(TARGET ${QSC_TESTS} PROPERTY CXX_STANDARD 11)
# Do not build the unit tests by default, since they take some time to build:
set_property(TARGET ${QSC_TESTS} PROPERTY EXCLUDE_FROM_ALL 1)



#enable_testing()
#
##add_test(unitTests unitTests)
#add_test(NAME unitTests COMMAND unitTests)
##add_test(mpiUnitTests tests/runMPIUnitTests)
#add_test(NAME mpiUnitTests COMMAND runMPIUnitTests WORKING_DIRECTORY tests)

# Set up "make test" so it runs the scripts I want it to run.
# See https://cmake.org/cmake/help/latest/command/add_custom_target.html
add_custom_target(test COMMAND ./runTests WORKING_DIRECTORY tests DEPENDS ${QSC_TESTS})

# Copy some files from the source to the build locations
message("CMAKE_CURRENT_SOURCE_DIR is ${CMAKE_CURRENT_SOURCE_DIR}")
message("CMAKE_CURRENT_BINARY_DIR is ${CMAKE_CURRENT_BINARY_DIR}")
set(TEST_FILES
  runTests
  runMPIUnitTests
  quasisymmetry_out.LandremanSengupta2019_section5.1.nc
  quasisymmetry_out.LandremanSengupta2019_section5.2.nc
  quasisymmetry_out.LandremanSengupta2019_section5.3.nc
  quasisymmetry_out.LandremanSengupta2019_section5.4.nc
  quasisymmetry_out.LandremanSengupta2019_section5.5.nc
  quasisymmetry_out.LandremanSenguptaPlunk_section5.1_order_r1_finite_r_linear.nc
  quasisymmetry_out.LandremanSenguptaPlunk_section5.2_order_r1_finite_r_nonlinear.nc
  quasisymmetry_out.LandremanSenguptaPlunk_section5.3_order_r1_finite_r_linear.nc
)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message("Building in root directory, so not copying files from source to build location.")
else()
  message("CMake version: ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")
  message("Copying files from source to build location.")
  foreach(TEST_FILE ${TEST_FILES})
    # It is preferable to create links, but CMake only introduced links in version 3.14
    if ((CMAKE_MAJOR_VERSION GREATER 3) OR (CMAKE_MINOR_VERSION GREATER 13))
      message("  linking ${TEST_FILE}")
      file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/tests/${TEST_FILE} ${CMAKE_CURRENT_BINARY_DIR}/tests/${TEST_FILE} COPY_ON_ERROR SYMBOLIC)
    else()
      message("  copying ${TEST_FILE}")
      file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests/${TEST_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests)
    endif()
  endforeach()
endif()


#message("GSL_FOUND = ${GSL_FOUND}")
#if(GSL_FOUND)
#  message("Found GSL")
#else()
#  message("Could not find GSL")
#endif()
