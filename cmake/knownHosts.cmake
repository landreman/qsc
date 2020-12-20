# Defaults for specific known hosts (i.e. HPC computing systems) can be put here. 

# For now, the following 4 lines are not used for anything. They might be useful in the future for detecting specific hosts.
cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
cmake_host_system_information(RESULT FQDN QUERY FQDN)
message("hostname is ${HOST}")
message("fully qualified domain name is ${FQDN}")


if(DEFINED ENV{TRAVIS_ARCH})
  message("Detected host is Travis-CI")
  # For Travis-CI, the flags "--mca btl_base_warn_component_unused 0 --mca orte_base_help_aggregate 0" are added to mpiexec
  # to avoid some warning messages.
  set(QSC_COMMAND_TO_SUBMIT_JOB "mpiexec -n NUM_PROCS --mca btl_base_warn_component_unused 0 --mca orte_base_help_aggregate 0")

elseif(DEFINED ENV{NERSC_HOST})
  message("Detected host is NERSC Cori")
  set(CMAKE_CXX_COMPILER "CC")

#elseif((DEFINED ENV{CLUSTER}) AND ($ENV{CLUSTER} STREQUAL "DRACO"))
elseif("$ENV{CLUSTER}" STREQUAL "DRACO")
  message("Detected host is IPP Draco")
  # CMake detects gnu compiler unless you specify Intel:
  set(CMAKE_CXX_COMPILER "mpiicpc")
  #set(BLA_VENDOR Intel10_64ilp_seq)

else()
  message("This host is not one with specific rules for qsc.")

endif()
