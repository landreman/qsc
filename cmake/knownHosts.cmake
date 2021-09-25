# Defaults for specific known hosts (i.e. HPC computing systems) can be put here. 

# For now, the following 4 lines are not used for anything. They might be useful in the future for detecting specific hosts.
cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
cmake_host_system_information(RESULT FQDN QUERY FQDN)
message("hostname is ${HOST}")
message("fully qualified domain name is ${FQDN}")


if(DEFINED ENV{GITHUB_ACTIONS})
  message("Detected host is Github Actions CI")
  # Some flags are added to mpiexec to avoid some warning messages.
  set(QSC_COMMAND_TO_SUBMIT_JOB "mpiexec -n NUM_PROCS --mca btl_base_warn_component_unused 0 --oversubscribe")
  # "--oversubscribe" is helpful because the macos environment is limited to 3 procs.
  #set(QSC_COMMAND_TO_SUBMIT_JOB "mpiexec -n NUM_PROCS --mca btl_base_warn_component_unused 0 --mca orte_base_help_aggregate 0")

elseif(DEFINED ENV{NERSC_HOST})
  message("Detected host is NERSC Cori")
  set(CMAKE_CXX_COMPILER "CC")

elseif("$ENV{CLUSTER}" STREQUAL "DRACO")
  message("Detected host is IPP Draco")
  # CMake detects gnu compiler unless you specify Intel:
  set(CMAKE_CXX_COMPILER "mpiicpc")
  #set(BLA_VENDOR Intel10_64ilp_seq)

elseif("$ENV{CLUSTER}" STREQUAL "COBRA")
  message("Detected host is IPP Cobra")
  # CMake detects gnu compiler unless you specify Intel:
  set(CMAKE_CXX_COMPILER "mpiicpc")
  # Compiled optimization flags recommended on the cobra website:
  # https://docs.mpcdf.mpg.de/doc/computing/software/compilers_languages.html#compiler-optimization-flags
  add_definitions(-xCORE-AVX512 -qopt-zmm-usage=high -ipo -O2)

elseif("$ENV{CLUSTER}" STREQUAL "RAVEN")
  message("Detected host is IPP Raven")
  # CMake detects gnu compiler unless you specify Intel:
  set(CMAKE_CXX_COMPILER "mpiicpc")

else()
  message("This host is not one with specific rules for qsc.")

endif()
