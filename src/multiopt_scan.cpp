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

  for (j = 0; j < N_FILTERS; j++) filters[j] = 0;
}

void MultiOptScan::scan() {
  big j_scan, jmod;
  int j, k, stage, stage_min, stage_max, index;
  std::valarray<int> indices;
  qscfloat val;

  indices.resize(ndim);
  for (j_scan = 0; j_scan < n_scan_all; j_scan++) {
    // From the global scan index j_scan, determine the indices for each dimension of the scan:
    jmod = j_scan;
    for (k = ndim - 1; k >= 0; k--) {
      indices[k] = jmod % params_n[k];
      jmod = (jmod - indices[k]) / params_n[k];
    }
    if (verbose > 0) {
      std::cout << "Scan indices:";
      for (j = 0; j < ndim; j++) std::cout << " " << indices[j];
      std::cout << std::endl;
    }

    // Initialize a MultiOpt with the default parameters:
    /*
    mo = MultiOpt();
    mo.opts.resize(mo_ref.opts.size());
    mo.opts[0].q = mo_ref.opts[0].q;
    for (j = 0; j < mo_ref.opts.size(); j++) {
      mo.opts[j].weight_iota = mo_ref.opts[j].weight_iota;
    }
    */
    mo = mo_ref;
    /*
    std::cout << "A mo_ref.opts[0].q.R0c[1]:" << mo_ref.opts[0].q.R0c[1] << std::endl;
    std::cout << "A     mo.opts[0].q.R0c[1]:" <<     mo.opts[0].q.R0c[1] << std::endl;
    mo_ref.opts[0].q.R0c[1] = 0.04;
    std::cout << "B mo_ref.opts[0].q.R0c[1]:" << mo_ref.opts[0].q.R0c[1] << std::endl;
    std::cout << "B     mo.opts[0].q.R0c[1]:" <<     mo.opts[0].q.R0c[1] << std::endl;
    */

    // Modify the MultiOpt parameters for this particular point in the scan:
    for (j = 0; j < ndim; j++) {
      val = params_vals[j][indices[j]];
      stage = params_stage[j];
      if (stage < 0) {
	stage_min = 0;
	stage_max = mo.opts.size();
      } else {
	stage_min = stage;
	stage_max = stage + 1;
      }
      
      for (stage = stage_min; stage < stage_max; stage++) {
	if (params[j].compare("R0c1") == 0) {
	  mo.opts[stage].q.R0c[1] = val;
	  if (verbose > 1) std::cout << "Setting R0c[1] for opt stage " << stage << " to " << val << std::endl;
	  
	} else if (params[j].compare("Z0s1") == 0) {
	  mo.opts[stage].q.Z0s[1] = val;
	  if (verbose > 1) std::cout << "Setting Z0s[1] for opt stage " << stage << " to " << val << std::endl;
	  
	} else if (params[j].compare("weight_iota") == 0) {
	  mo.opts[stage].weight_iota = val;
	  if (verbose > 1) std::cout << "Setting weight_iota for opt stage " << stage << " to " << val << std::endl;
	  
	} else if (params[j].compare("weight_B20") == 0) {
	  mo.opts[stage].weight_B20 = val;
	  if (verbose > 1) std::cout << "Setting weight_B20 for opt stage " << stage << " to " << val << std::endl;
	  
	} else if (params[j].compare("weight_r_singularity") == 0) {
	  mo.opts[stage].weight_r_singularity = val;
	  if (verbose > 1) std::cout << "Setting weight_r_singularity for opt stage " << stage << " to " << val << std::endl;
	  
	} else if (params[j].compare("weight_grad_grad_B") == 0) {
	  mo.opts[stage].weight_grad_grad_B = val;
	  if (verbose > 1) std::cout << "Setting weight_grad_grad_B for opt stage " << stage << " to " << val << std::endl;
	  
	} else if (params[j].compare("weight_axis_length") == 0) {
	  mo.opts[stage].weight_axis_length = val;
	  if (verbose > 1) std::cout << "Setting weight_axis_length for opt stage " << stage << " to " << val << std::endl;
	  
	} else if (params[j].compare("weight_standard_deviation_of_R") == 0) {
	  mo.opts[stage].weight_standard_deviation_of_R = val;
	  if (verbose > 1) std::cout << "Setting weight_standard_deviation_of_R for opt stage " << stage << " to " << val << std::endl;
	  
	} else {
	  std::cout << "Unrecognized entry in params: " << params[j] << std::endl;
	  throw std::runtime_error("Unrecognized entry in params");
	}
      }
    }

    // Run the MultiOpt:
    mo.optimize();

    // See if the final configuration passes the filters:
    if (!keep_all) {
      index = mo.opts.size() - 1;
      
      if (mo.opts[index].q.grid_min_R0 < min_R0_to_keep) {
	filters[REJECTED_DUE_TO_R0]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to min_R0." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed min_R0 filter." << std::endl;
      
      if (std::abs(mo.opts[index].q.iota) < min_iota_to_keep) {
	filters[REJECTED_DUE_TO_IOTA]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to iota." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed iota filter." << std::endl;
      
      if (mo.opts[index].q.grid_max_elongation > max_elongation_to_keep) {
	filters[REJECTED_DUE_TO_ELONGATION]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to elongation." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed elongation filter." << std::endl;
      
      if (mo.opts[index].q.grid_min_L_grad_B < min_L_grad_B_to_keep) {
	filters[REJECTED_DUE_TO_L_GRAD_B]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to L_grad_B." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed L_grad_B filter." << std::endl;
      
      if (mo.opts[index].q.B20_grid_variation > max_B20_variation_to_keep) {
	filters[REJECTED_DUE_TO_B20_VARIATION]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to B20 variation." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed B20 variation filter." << std::endl;
      
      if (mo.opts[index].q.d2_volume_d_psi2 > max_d2_volume_d_psi2_to_keep) {
	filters[REJECTED_DUE_TO_D2_VOLUME_D_PSI2]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to d2_volume_d_psi2." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed d2_volume_d_psi2 filter." << std::endl;
      
      if (mo.opts[index].q.DMerc_times_r2 < min_DMerc_to_keep) {
	filters[REJECTED_DUE_TO_DMERC]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to DMerc." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed DMerc filter." << std::endl;
      
      if (mo.opts[index].q.grid_min_L_grad_grad_B < min_L_grad_grad_B_to_keep) {
	filters[REJECTED_DUE_TO_L_GRAD_GRAD_B]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to L_grad_grad_B." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed L_grad_grad_B filter." << std::endl;
      
      if (mo.opts[index].q.r_singularity_robust < min_r_singularity_to_keep) {
	filters[REJECTED_DUE_TO_R_SINGULARITY]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to r_singularity." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed r_singularity filter." << std::endl;
      
    }
  }
}
