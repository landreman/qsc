#include <mpi.h>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iomanip>
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
  write_netcdf();
}

void MultiOptScan::init() {
  int j, k;
  
  // Set the values for each dimension of the scan.
  ndim = params.size();
  assert (params_max.size() == ndim);
  assert (params_min.size() == ndim);
  assert (params_n.size() == ndim);
  assert (params_log.size() == ndim);
  assert (params_stage.size() == ndim);
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
  axis_nmax_plus_1 = mo_ref.opts[0].q.R0c.size();
  for (j = 0; j < mo_ref.opts.size(); j++) axis_nmax_plus_1 += mo_ref.opts[j].fourier_refine;

}

void MultiOptScan::scan() {
  big j_scan, jmod;
  int j, k, stage, stage_min, stage_max, index;
  std::valarray<int> indices;
  qscfloat val;

  scan_eta_bar.resize(n_scan_all, 0.0);
  scan_sigma0.resize(n_scan_all, 0.0);
  scan_B2c.resize(n_scan_all, 0.0);
  scan_B2s.resize(n_scan_all, 0.0);
  scan_min_R0.resize(n_scan_all, 0.0);
  scan_max_curvature.resize(n_scan_all, 0.0);
  scan_iota.resize(n_scan_all, 0.0);
  scan_max_elongation.resize(n_scan_all, 0.0);
  scan_min_L_grad_B.resize(n_scan_all, 0.0);
  scan_min_L_grad_grad_B.resize(n_scan_all, 0.0);
  scan_r_singularity.resize(n_scan_all, 0.0);
  scan_d2_volume_d_psi2.resize(n_scan_all, 0.0);
  scan_DMerc_times_r2.resize(n_scan_all, 0.0);
  scan_B20_variation.resize(n_scan_all, 0.0);
  scan_B20_residual.resize(n_scan_all, 0.0);
  scan_standard_deviation_of_R.resize(n_scan_all, 0.0);
  scan_standard_deviation_of_Z.resize(n_scan_all, 0.0);

  scan_helicity.resize(n_scan_all, 0);

  scan_R0c.resize(axis_nmax_plus_1, n_scan_all, 0.0);
  scan_R0s.resize(axis_nmax_plus_1, n_scan_all, 0.0);
  scan_Z0c.resize(axis_nmax_plus_1, n_scan_all, 0.0);
  scan_Z0s.resize(axis_nmax_plus_1, n_scan_all, 0.0);

  scan_weight_B20.resize(n_scan_all, 0.0);
  scan_weight_iota.resize(n_scan_all, 0.0);
  scan_target_iota.resize(n_scan_all, 0.0);
  scan_weight_elongation.resize(n_scan_all, 0.0);
  scan_weight_curvature.resize(n_scan_all, 0.0);
  scan_weight_R0.resize(n_scan_all, 0.0);
  scan_target_min_R0.resize(n_scan_all, 0.0);
  scan_weight_d2_volume_d_psi2.resize(n_scan_all, 0.0);
  scan_max_d2_volume_d_psi2.resize(n_scan_all, 0.0);
  scan_weight_XY2.resize(n_scan_all, 0.0);
  scan_weight_XY2Prime.resize(n_scan_all, 0.0);
  scan_weight_Z2.resize(n_scan_all, 0.0);
  scan_weight_Z2Prime.resize(n_scan_all, 0.0);
  scan_weight_XY3.resize(n_scan_all, 0.0);
  scan_weight_XY3Prime.resize(n_scan_all, 0.0);
  scan_weight_grad_B.resize(n_scan_all, 0.0);
  scan_weight_grad_grad_B.resize(n_scan_all, 0.0);
  scan_weight_r_singularity.resize(n_scan_all, 0.0);
  scan_weight_axis_length.resize(n_scan_all, 0.0);
  scan_weight_standard_deviation_of_R.resize(n_scan_all, 0.0);
  
  indices.resize(ndim);
  n_scan = 0;
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

    // Index of the last optimization stage in mo.opts:
    index = mo.opts.size() - 1;

    // See if the final configuration passes the filters:
    if (!keep_all) {
      
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

    // If we make it this far, then save the configuration.
    scan_eta_bar[n_scan] = mo.opts[index].q.eta_bar;
    scan_sigma0[n_scan] = mo.opts[index].q.sigma0;
    scan_B2c[n_scan] = mo.opts[index].q.B2c;
    scan_B2s[n_scan] = mo.opts[index].q.B2s;
    scan_min_R0[n_scan] = mo.opts[index].q.grid_min_R0;
    scan_max_curvature[n_scan] = mo.opts[index].q.grid_max_curvature;
    scan_iota[n_scan] = mo.opts[index].q.iota;
    scan_max_elongation[n_scan] = mo.opts[index].q.grid_max_elongation;
    scan_min_L_grad_B[n_scan] = mo.opts[index].q.grid_min_L_grad_B;
    scan_min_L_grad_grad_B[n_scan] = mo.opts[index].q.grid_min_L_grad_grad_B;
    scan_r_singularity[n_scan] = mo.opts[index].q.r_singularity_robust;
    scan_d2_volume_d_psi2[n_scan] = mo.opts[index].q.d2_volume_d_psi2;
    scan_DMerc_times_r2[n_scan] = mo.opts[index].q.DMerc_times_r2;
    scan_B20_variation[n_scan] = mo.opts[index].q.B20_grid_variation;
    scan_B20_residual[n_scan] = mo.opts[index].q.B20_residual;
    scan_standard_deviation_of_R[n_scan] = mo.opts[index].q.standard_deviation_of_R;
    scan_standard_deviation_of_Z[n_scan] = mo.opts[index].q.standard_deviation_of_Z;
    scan_helicity[n_scan] = mo.opts[index].q.helicity;
    for (j = 0; j < axis_nmax_plus_1; j++) {
      scan_R0c(j, n_scan) = mo.opts[index].q.R0c[j];
      scan_R0s(j, n_scan) = mo.opts[index].q.R0s[j];
      scan_Z0c(j, n_scan) = mo.opts[index].q.Z0c[j];
      scan_Z0s(j, n_scan) = mo.opts[index].q.Z0s[j];
    }

    scan_weight_B20[n_scan] = mo.opts[index].weight_B20;
    scan_weight_iota[n_scan] = mo.opts[index].weight_iota;
    scan_target_iota[n_scan] = mo.opts[index].target_iota;
    scan_weight_elongation[n_scan] = mo.opts[index].weight_elongation;
    scan_weight_curvature[n_scan] = mo.opts[index].weight_curvature;
    scan_weight_R0[n_scan] = mo.opts[index].weight_R0;
    scan_target_min_R0[n_scan] = mo.opts[index].min_R0;
    scan_weight_d2_volume_d_psi2[n_scan] = mo.opts[index].weight_d2_volume_d_psi2;
    scan_max_d2_volume_d_psi2[n_scan] = mo.opts[index].max_d2_volume_d_psi2;
    scan_weight_XY2[n_scan] = mo.opts[index].weight_XY2;
    scan_weight_XY2Prime[n_scan] = mo.opts[index].weight_XY2Prime;
    scan_weight_Z2[n_scan] = mo.opts[index].weight_Z2;
    scan_weight_Z2Prime[n_scan] = mo.opts[index].weight_Z2Prime;
    scan_weight_XY3[n_scan] = mo.opts[index].weight_XY3;
    scan_weight_XY3Prime[n_scan] = mo.opts[index].weight_XY3Prime;
    scan_weight_grad_B[n_scan] = mo.opts[index].weight_grad_B;
    scan_weight_grad_grad_B[n_scan] = mo.opts[index].weight_grad_grad_B;
    scan_weight_r_singularity[n_scan] = mo.opts[index].weight_r_singularity;
    scan_weight_axis_length[n_scan] = mo.opts[index].weight_axis_length;
    scan_weight_standard_deviation_of_R[n_scan] = mo.opts[index].weight_standard_deviation_of_R;
    
    n_scan += 1;
  }

  filters[ATTEMPTS] = n_scan_all;
  filters[KEPT] = n_scan;
  big total_rejected = 0;
  for (j = REJECTED_DUE_TO_R0; j < N_FILTERS; j++) total_rejected += filters[j];
  for (j = 0; j < N_FILTERS; j++) filter_fractions[j] = ((qscfloat)filters[j]) / filters[0];

  int width = 13;
  std::cout << std::setprecision(4) << std::endl;
  std::cout << "Summary of scan results:                           (fractions in parentheses)" << std::endl;
  std::cout << "  Configurations attempted:          " << std::setw(width) << filters[ATTEMPTS] << std::endl;
  std::cout << "  Rejected due to min_R0:            " << std::setw(width) << filters[REJECTED_DUE_TO_R0]
	    << " (" << filter_fractions[REJECTED_DUE_TO_R0] << ")" << std::endl;
  std::cout << "  Rejected due to min iota:          " << std::setw(width) << filters[REJECTED_DUE_TO_IOTA]
	    << " (" << filter_fractions[REJECTED_DUE_TO_IOTA] << ")" << std::endl;
  std::cout << "  Rejected due to max elongation:    " << std::setw(width) << filters[REJECTED_DUE_TO_ELONGATION]
	    << " (" << filter_fractions[REJECTED_DUE_TO_ELONGATION] << ")" << std::endl;
  std::cout << "  Rejected due to min L_grad_B:      " << std::setw(width) << filters[REJECTED_DUE_TO_L_GRAD_B]
	    << " (" << filter_fractions[REJECTED_DUE_TO_L_GRAD_B] << ")" << std::endl;
  std::cout << "  Rejected due to B20 variation:     " << std::setw(width) << filters[REJECTED_DUE_TO_B20_VARIATION]
	    << " (" << filter_fractions[REJECTED_DUE_TO_B20_VARIATION] << ")" << std::endl;
  std::cout << "  Rejected due to min L_grad_grad_B: " << std::setw(width) << filters[REJECTED_DUE_TO_L_GRAD_GRAD_B]
	    << " (" << filter_fractions[REJECTED_DUE_TO_L_GRAD_GRAD_B] << ")" << std::endl;
  std::cout << "  Rejected due to d2_volume_d_psi2:  " << std::setw(width) << filters[REJECTED_DUE_TO_D2_VOLUME_D_PSI2]
	    << " (" << filter_fractions[REJECTED_DUE_TO_D2_VOLUME_D_PSI2] << ")" << std::endl;
  std::cout << "  Rejected due to DMerc:             " << std::setw(width) << filters[REJECTED_DUE_TO_DMERC]
	    << " (" << filter_fractions[REJECTED_DUE_TO_DMERC] << ")" << std::endl;
  std::cout << "  Rejected due to r_singularity:     " << std::setw(width) << filters[REJECTED_DUE_TO_R_SINGULARITY]
	    << " (" << filter_fractions[REJECTED_DUE_TO_R_SINGULARITY] << ")" << std::endl;
  std::cout << "  Total rejected:                    " << std::setw(width) << total_rejected
	    << " (" << ((qscfloat)total_rejected) / filters[ATTEMPTS] << ")" << std::endl;
  std::cout << "  Kept:                              " << std::setw(width) << n_scan
	    << " (" << filter_fractions[KEPT] << ")" << std::endl;
  std::cout << "  Kept + rejected:                   " << std::setw(width) << n_scan + total_rejected
	    << std::endl;
  
  /*
    std::cout << std::setprecision(4) << "Time elapsed, summed over processors:" << std::endl;
    width = 10;
    
    std::cout << "  Time for random number generation: " << std::setw(width) << timing[TIME_RANDOM]
    << " (" << timing[TIME_RANDOM] / timing_total << ")" << std::endl;
    
    std::cout << "  Time for init_axis:                " << std::setw(width) << timing[TIME_INIT_AXIS]
    << " (" << timing[TIME_INIT_AXIS] / timing_total << ")" << std::endl;
    
    std::cout << "  Time for solving sigma equation:   " << std::setw(width) << timing[TIME_SIGMA_EQUATION]
    << " (" << timing[TIME_SIGMA_EQUATION] / timing_total << ")" << std::endl;
    
    std::cout << "  Time for O(r^1) diagnostics:       " << std::setw(width) << timing[TIME_R1_DIAGNOSTICS]
    << " (" << timing[TIME_R1_DIAGNOSTICS] / timing_total << ")" << std::endl;
    
    if (q.at_least_order_r2) {
      std::cout << "  Time for calculate_r2:             " << std::setw(width) << timing[TIME_CALCULATE_R2]
		<< " (" << timing[TIME_CALCULATE_R2] / timing_total << ")" << std::endl;
    
      std::cout << "  Time for mercier:                  " << std::setw(width) << timing[TIME_MERCIER]
		<< " (" << timing[TIME_MERCIER] / timing_total << ")" << std::endl;
    
      std::cout << "  Time for grad grad B tensor:       " << std::setw(width) << timing[TIME_GRAD_GRAD_B_TENSOR]
		<< " (" << timing[TIME_GRAD_GRAD_B_TENSOR] / timing_total << ")" << std::endl;
    
      std::cout << "  Time for r_singularity:            " << std::setw(width) << timing[TIME_R_SINGULARITY]
		<< " (" << timing[TIME_R_SINGULARITY] / timing_total << ")" << std::endl;
    }
  */
  if (n_scan < 1000) {
    std::cout << std::setprecision(2) << std::endl;
    std::cout << "min_R0: " << scan_min_R0 << std::endl;
    std::cout << std::endl;
    std::cout << "iota: " << scan_iota << std::endl;
    std::cout << std::endl;
    std::cout << "max_elongation: " << scan_max_elongation << std::endl;
    std::cout << std::endl;
    std::cout << "helicity:";
    for (j = 0; j < n_scan; j++) std::cout << " " << scan_helicity[j];
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "min_L_grad_B: " << scan_min_L_grad_B << std::endl;
    if (mo_ref.opts[0].q.at_least_order_r2) {
      std::cout << std::endl;
      std::cout << "min_L_grad_grad_B: " << scan_min_L_grad_grad_B << std::endl;
      std::cout << std::endl;
      std::cout << "d2_volume_d_psi2: " << scan_d2_volume_d_psi2 << std::endl;
      std::cout << std::endl;
      std::cout << "r_singularity: " << scan_r_singularity << std::endl;
      std::cout << std::endl;
      std::cout << "B20_variation: " << scan_B20_variation << std::endl;
    }
    std::cout << std::endl;
    
    // Restore precision for printing
    if (single) {
      std::cout.precision(8);
      std::cerr.precision(8);
    } else {
      std::cout.precision(15);
      std::cerr.precision(15);
    }
  }
}
