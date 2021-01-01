#include <mpi.h>
#include "scan.hpp"

using namespace qsc;

void Scan::defaults() {
  // Set defaults.
  mpi_comm = MPI_COMM_WORLD;
  max_seconds = 60;
  max_keep_per_proc = 1000;
  max_attempts_per_proc = 10000;
  deterministic = false;
  
  eta_bar_min = 1.0;
  eta_bar_max = 1.0;
  sigma0_min = 0.0;
  sigma0_max = 0.0;
  B2c_min = 0.0;
  B2c_max = 0.0;
  B2s_min = 0.0;
  B2s_max = 0.0;
    
  R0c_min.resize(1, 1.0);
  R0c_max.resize(1, 1.0);
  R0s_min.resize(1, 0.0);
  R0s_max.resize(1, 0.0);
  Z0c_min.resize(1, 0.0);
  Z0c_max.resize(1, 0.0);
  Z0s_min.resize(1, 0.0);
  Z0s_max.resize(1, 0.0);

  eta_bar_scan_option = SCAN_OPTION_LINEAR;
  sigma0_scan_option = SCAN_OPTION_LINEAR;
  B2c_scan_option = SCAN_OPTION_LINEAR;
  B2s_scan_option = SCAN_OPTION_LINEAR;
  fourier_scan_option = SCAN_OPTION_LINEAR;

  keep_all = true;
  min_R0_to_keep = -1.0;
  min_iota_to_keep = -1.0;
  max_elongation_to_keep = 10.0;
  min_L_grad_B_to_keep = -1.0;
  min_L_grad_grad_B_to_keep = -1.0;
  max_B20_variation_to_keep = 1.0e+30;
  min_r_singularity_to_keep = -1.0;
  max_d2_volume_d_psi2_to_keep = 1.0e+30;
  min_DMerc_to_keep = -1.0e+30;
}

Scan::Scan() {
  defaults();
}

void Scan::run(std::string directory_and_infile) {
  std::string directory_and_outfile = qsc::outfile(directory_and_infile);
  input(directory_and_infile);
  random();
  // calculate();
  // write_netcdf(directory_and_outfile);
}
