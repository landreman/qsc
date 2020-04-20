# Defaults for specific known hosts (i.e. HPC computing systems) can be put here. 

# For now, the following 4 lines are not used for anything. They might be useful in the future for detecting specific hosts.
cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
cmake_host_system_information(RESULT FQDN QUERY FQDN)
message("hostname is ${HOST}")
message("fully qualified domain name is ${FQDN}")


if(DEFINED ENV{TRAVIS_ARCH})
  message("Detected host is Travis-CI")
  set(QSC_COMMAND_TO_SUBMIT_JOB "mpiexec -n NUM_PROCS --mca btl_base_warn_component_unused 0 --mca orte_base_help_aggregate 0")

elseif(DEFINED ENV{NERSC_HOST})
  message("Detected host is NERSC Cori")

else()
  message("This host is not one already known to qsc.")

endif()
