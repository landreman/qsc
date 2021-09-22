#include <mpi.h>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include "multiopt_scan.hpp"

using namespace qsc;

void MultiOptScan::defaults() {
  // Set defaults.
  mpi_comm = MPI_COMM_WORLD;
  verbose = 1;
  max_seconds = 60;
  save_period = 60;
  max_keep_per_proc = 1000;
  max_attempts_per_proc = -1;
  
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

MultiOptScan::MultiOptScan() {
  defaults();
}

void MultiOptScan::run(std::string directory_and_infile) {
  outfilename = qsc::outfile(directory_and_infile);
  input(directory_and_infile);
  init();
  scan();
  //write_netcdf();
}

void MultiOptScan::init() {
  int j, k;
  
  // Set the values for each dimension of the scan.
  ndim = params.size();
  assert (params_max.size() == ndim);
  assert (params_min.size() == ndim);
  assert (params_n.size() == ndim);
  assert (params_log.size() == ndim);
  params_vals.resize(ndim);
  n_scan_all = 1;
  for (j = 0; j < ndim; j++) {
    if (params_n[j] < 2) throw std::runtime_error("params_n values must all be >= 2");
    n_scan_all *= params_n[j];
    params_vals[j].resize(params_n[j]);
    if (params_log[j]) {
      // Logarithmic spacing
      for (int k = 0; k < params_n[j]; k++)
	params_vals[j][k] = exp(log(params_min[j]) + (log(params_max[j]) - log(params_min[j])) * k / (params_n[j] - 1));
    } else {
      // Linear spacing
      for (int k = 0; k < params_n[j]; k++)
	params_vals[j][k] = params_min[j] + (params_max[j] - params_min[j]) * k / (params_n[j] - 1);
    }
    if (verbose > 0)
      std::cout << "Values for parameter " << params[j] << ": " << params_vals[j] << std::endl;
  }
  if (verbose > 0) std::cout << "Total number of points in scan: " << n_scan_all << std::endl;
}

void MultiOptScan::scan() {
}
