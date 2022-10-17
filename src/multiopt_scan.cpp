#include <mpi.h>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include "multiopt_scan.hpp"

using namespace qsc;

void MultiOptScan::defaults() {
  // Set defaults.
  mpi_comm = MPI_COMM_WORLD;
  verbose = 1;
  max_seconds = 60;
  print_status_period = 20.0;
  save_period = 5 * 60.0;
  quit_after_init = false;
  
  keep_all = true;
  min_R0_to_keep = -1.0;
  min_iota_to_keep = -1.0;
  max_elongation_to_keep = 10.0;
  min_L_grad_B_to_keep = -1.0;
  min_L_grad_grad_B_to_keep = -1.0;
  max_B20_variation_to_keep = 1.0e+30;
  min_r_singularity_to_keep = -1.0;
  max_d2_volume_d_psi2_to_keep = 1.0e+30;
  min_DMerc_times_r2_to_keep = -1.0e+30;
  max_XY2_to_keep = 1.0e+30;
  max_Z2_to_keep = 1.0e+30;
  max_XY3_to_keep = 1.0e+30;
  max_d_XY2_d_varphi_to_keep = 1.0e+30;
  max_d_Z2_d_varphi_to_keep = 1.0e+30;
  max_d_XY3_d_varphi_to_keep = 1.0e+30;
}

MultiOptScan::MultiOptScan() {
  defaults();
}

void MultiOptScan::run(std::string directory_and_infile) {
  outfilename = qsc::outfile(directory_and_infile);
  input(directory_and_infile);
  init();
  if (!quit_after_init) {
    scan();
    write_netcdf();
  }
}

void MultiOptScan::init() {
  int j, k;
  
  // Initialize MPI-related variables
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Comm_size(mpi_comm, &n_procs);
  proc0 = (mpi_rank == 0);
  MPI_Barrier(mpi_comm);

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
    if (verbose > -1 && proc0)
      std::cout << "Values for parameter " << params[j] << ": " << params_vals[j] << std::endl;
  }
  if (verbose > -1 && proc0) std::cout << "Total number of points in scan: " << n_scan_all << std::endl;

  axis_nmax_plus_1 = mo_ref.opts[0].q.R0c.size();
  for (j = 0; j < mo_ref.opts.size(); j++) axis_nmax_plus_1 += mo_ref.opts[j].fourier_refine;
  int n_fourier_parameters = axis_nmax_plus_1 * 4;

  n_parameters = n_parameters_base + n_fourier_parameters;
  parameters_single.resize(n_parameters);
  int_parameters_single.resize(n_int_parameters);
  if (proc0) {
    n_solves_kept.resize(n_procs, 0);
    attempts_per_proc.resize(n_procs, 0);
    for (j = 0; j < N_FILTERS; j++) filters[j] = 0;
  }
  total_cpu_seconds = 0.0;
  n_evals = 0;

}

void MultiOptScan::scan() {
  MPI_Status mpi_status;
  int intbuf[1];
  int proc_that_finished, j;
  std::chrono::time_point<std::chrono::steady_clock> start_time2, end_time, start_time_print_status, start_time_save;
  std::chrono::duration<double> elapsed;
  
  if (n_procs < 2)
    throw std::runtime_error("For MultiOptScan, the number of MPI processes must be at least 2.");
  if (n_procs - 1 > n_scan_all)
    throw std::runtime_error("For MultiOptScan, the number of MPI processes cannot exceed n_scan_all + 1.");
  
  start_time = std::chrono::steady_clock::now();
  start_time_print_status = start_time;
  start_time_save = start_time;
  if (proc0) {
    // Allocate the big arrays for global results:
    if (verbose > 0) std::cout << "proc0 about to allocate big arrays." << std::endl;
    start_time2 = std::chrono::steady_clock::now();
    parameters.resize(n_parameters, n_scan_all, 0.0);
    int_parameters.resize(n_int_parameters_base * n_scan_all, 0);
    end_time = std::chrono::steady_clock::now();
    elapsed = end_time - start_time2;
    std::cout << "Time for big allocations: " << elapsed.count() << " seconds" << std::endl;
    
    // Send off the first batch of work:
    for (j = 1; j < n_procs; j++) {
      intbuf[0] = j - 1;
      MPI_Send(intbuf, 1, MPI_INT, j, 0, mpi_comm);
    }
    // Handle the rest of the tasks:
    for (j = n_procs - 1; j < n_scan_all; j++) {
      // Print a status summary every time interval print_status_period;
      end_time = std::chrono::steady_clock::now();
      elapsed = end_time - start_time_print_status;
      if (elapsed.count() > print_status_period) {
	print_status();
	start_time_print_status = end_time;
      }
      
      // Save intermediate results every time interval save_period;
      end_time = std::chrono::steady_clock::now();
      elapsed = end_time - start_time_save;
      if (elapsed.count() > save_period) {
	print_status();
	filter_global_arrays();
	write_netcdf();
	start_time_save = end_time;
      }
	
      // Wait for any proc to finish:
      proc_that_finished = proc0_recv();
      //MPI_Recv(intbuf, 1, MPI_INT, MPI_ANY_SOURCE, 0, mpi_comm, &mpi_status);
      //proc_that_finished = intbuf[0];
      if (verbose > 0) std::cout << "proc0: I see that proc " << proc_that_finished << " just finished, so I'll send it index " << j << std::endl;
      intbuf[0] = j;
      // Give that proc the next piece of work:
      MPI_Send(intbuf, 1, MPI_INT, proc_that_finished, 0, mpi_comm);
    }
    // Receive the final completion messages from all other procs:
    for (j = 1; j < n_procs; j++) {
      //MPI_Recv(intbuf, 1, MPI_INT, MPI_ANY_SOURCE, 0, mpi_comm, &mpi_status);
      //proc_that_finished = intbuf[0];
      proc_that_finished = proc0_recv();
      if (verbose > 0) std::cout << "proc0: I see that proc " << proc_that_finished << " just finished. No further work." << std::endl;
      // Tell that proc to exit:
      intbuf[0] = -1;
      MPI_Send(intbuf, 1, MPI_INT, j, 0, mpi_comm);
    }
    
  } else {
    // Code for procs >0.

    bool keep_going = true;
    while (keep_going) {
      // Wait for an instruction from proc0:
      MPI_Recv(intbuf, 1, MPI_INT, 0, 0, mpi_comm, &mpi_status);
      if (intbuf[0] < 0) {
	keep_going = false;
	if (verbose > 0) std::cout << "proc " << mpi_rank << " is exiting." << std::endl;
      } else {
	if (verbose > 0) std::cout << "proc " << mpi_rank << " is handling index " << intbuf[0] << std::endl;
	eval_scan_index(intbuf[0]);
	// Return results to proc0.
	if (verbose > 0) std::cout << "proc " << mpi_rank << " is returning results for index " << intbuf[0] << std::endl;
	MPI_Send(&int_parameters_single[0], n_int_parameters,      MPI_INT, 0, 0, mpi_comm);
	MPI_Send(    &parameters_single[0],     n_parameters, MPI_QSCFLOAT, 0, 0, mpi_comm);
	//intbuf[0] = mpi_rank;
	//MPI_Send(intbuf, 1, MPI_INT, 0, 0, mpi_comm);
      }
    }
  }
  
  end_time = std::chrono::steady_clock::now();
  elapsed = end_time - start_time;
  std::cout << "Proc " << mpi_rank << " finished after " << elapsed.count() << " seconds" << std::endl;

  MPI_Barrier(mpi_comm);
  if (proc0) {
    print_status();
    filter_global_arrays();
  }
}

/**
 * In this routine, proc0 receives data from another proc, and stores it in the global arrays.
 * This routine is only ever run on proc0.
 */
int MultiOptScan::proc0_recv() {
  MPI_Status mpi_status;
  int proc_that_finished, j, j_scan;
  
  MPI_Recv(&int_parameters_single[0], n_int_parameters,      MPI_INT,     MPI_ANY_SOURCE, 0, mpi_comm, &mpi_status);
  proc_that_finished = mpi_status.MPI_SOURCE;
  MPI_Recv(    &parameters_single[0],     n_parameters, MPI_QSCFLOAT, proc_that_finished, 0, mpi_comm, &mpi_status);

  // Store results in the global arrays:
  j_scan = int_parameters_single[0];
  for (j = 0; j < n_parameters; j++)
    parameters(j, j_scan) = parameters_single[j];
  for (j = 0; j < n_int_parameters_base; j++)
    int_parameters[j + n_int_parameters_base * j_scan] = int_parameters_single[j];
  for (j = 0; j < N_FILTERS; j++) {
    filters_local[j] = int_parameters_single[j + n_int_parameters_base];
    filters[j] += filters_local[j];
  }
  attempts_per_proc[proc_that_finished]++;
  if (filters_local[KEPT] > 0) n_solves_kept[proc_that_finished]++;
  total_cpu_seconds += parameters_single[0];
  n_evals += int_parameters_single[2];
  
  return proc_that_finished;
}

/**
 * Evaluate a single point in the parameter scan.
 * This routine is only run on procs > 0.
 */
void MultiOptScan::eval_scan_index(int j_scan) {
  big jmod;
  int j, k, stage, stage_min, stage_max, index;
  std::valarray<int> indices;
  qscfloat val;
  std::chrono::time_point<std::chrono::steady_clock> start_time_single, end_time_single;
  std::chrono::duration<double> elapsed;
  bool passed_filters;
  
  // Some indices of interest:
  // j_scan is the index into the complete tensor-product set of points.
  // filters_local[ATTEMPTS] is the number of configurations considered so far on this process.
  
  indices.resize(ndim);
  start_time_single = std::chrono::steady_clock::now();
  for (j = 0; j < N_FILTERS; j++) filters_local[j] = 0;
  filters_local[ATTEMPTS]++;

  // From the global scan index j_scan, determine the indices for each dimension of the scan:
  jmod = j_scan;
  for (k = ndim - 1; k >= 0; k--) {
    indices[k] = jmod % params_n[k];
    jmod = (jmod - indices[k]) / params_n[k];
  }
  if (verbose > 0 && proc0) {
    std::cout << "Scan indices:";
    for (j = 0; j < ndim; j++) std::cout << " " << indices[j];
    std::cout << std::endl;
  }

  // Initialize a MultiOpt with the default parameters:
  mo = mo_ref;
  
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
	
      } else if (params[j].compare("target_iota") == 0) {
	mo.opts[stage].target_iota = val;
	if (verbose > 1) std::cout << "Setting target_iota for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_elongation") == 0) {
	mo.opts[stage].weight_elongation = val;
	if (verbose > 1) std::cout << "Setting weight_elongation for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_B20") == 0) {
	mo.opts[stage].weight_B20 = val;
	if (verbose > 1) std::cout << "Setting weight_B20 for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_r_singularity") == 0) {
	mo.opts[stage].weight_r_singularity = val;
	if (verbose > 1) std::cout << "Setting weight_r_singularity for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_grad_B") == 0) {
	mo.opts[stage].weight_grad_B = val;
	if (verbose > 1) std::cout << "Setting weight_grad_B for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_grad_grad_B") == 0) {
	mo.opts[stage].weight_grad_grad_B = val;
	if (verbose > 1) std::cout << "Setting weight_grad_grad_B for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_axis_length") == 0) {
	mo.opts[stage].weight_axis_length = val;
	if (verbose > 1) std::cout << "Setting weight_axis_length for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("target_axis_length") == 0) {
	mo.opts[stage].target_axis_length = val;
	if (verbose > 1) std::cout << "Setting target_axis_length for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_standard_deviation_of_R") == 0) {
	mo.opts[stage].weight_standard_deviation_of_R = val;
	if (verbose > 1) std::cout << "Setting weight_standard_deviation_of_R for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_XY2") == 0) {
	mo.opts[stage].weight_XY2 = val;
	if (verbose > 1) std::cout << "Setting weight_XY2 for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_Z2") == 0) {
	mo.opts[stage].weight_Z2 = val;
	if (verbose > 1) std::cout << "Setting weight_Z2 for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_XY3") == 0) {
	mo.opts[stage].weight_XY3 = val;
	if (verbose > 1) std::cout << "Setting weight_XY3 for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_XY2Prime") == 0) {
	mo.opts[stage].weight_XY2Prime = val;
	if (verbose > 1) std::cout << "Setting weight_XY2Prime for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_XY2PrimePrime") == 0) {
	mo.opts[stage].weight_XY2PrimePrime = val;
	if (verbose > 1) std::cout << "Setting weight_XY2PrimePrime for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_Z2Prime") == 0) {
	mo.opts[stage].weight_Z2Prime = val;
	if (verbose > 1) std::cout << "Setting weight_Z2Prime for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_XY3Prime") == 0) {
	mo.opts[stage].weight_XY3Prime = val;
	if (verbose > 1) std::cout << "Setting weight_XY3Prime for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_XY3PrimePrime") == 0) {
	mo.opts[stage].weight_XY3PrimePrime = val;
	if (verbose > 1) std::cout << "Setting weight_XY3PrimePrime for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_B20_mean") == 0) {
	mo.opts[stage].weight_B20_mean = val;
	if (verbose > 1) std::cout << "Setting weight_B20_mean for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_d2_volume_d_psi2") == 0) {
	mo.opts[stage].weight_d2_volume_d_psi2 = val;
	if (verbose > 1) std::cout << "Setting weight_d2_volume_d_psi2 for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("max_d2_volume_d_psi2") == 0) {
	mo.opts[stage].max_d2_volume_d_psi2 = val;
	if (verbose > 1) std::cout << "Setting max_d2_volume_d_psi2 for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("weight_DMerc_times_r2") == 0) {
	mo.opts[stage].weight_DMerc_times_r2 = val;
	if (verbose > 1) std::cout << "Setting weight_DMerc_times_r2 for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("min_DMerc_times_r2") == 0) {
	mo.opts[stage].min_DMerc_times_r2 = val;
	if (verbose > 1) std::cout << "Setting min_DMerc_times_r2 for opt stage " << stage << " to " << val << std::endl;
	
      } else if (params[j].compare("eta_bar") == 0) {
	mo.opts[stage].q.eta_bar = val;
	if (verbose > 1) std::cout << "Setting eta_bar for opt stage " << stage << " to " << val << std::endl;
	if (stage != 0) throw std::runtime_error("In a multiopt_scan, eta_bar should only be set at stage 0.");
	
      } else if (params[j].compare("p2") == 0) {
	mo.opts[stage].q.p2 = val;
	if (verbose > 1) std::cout << "Setting p2 for opt stage " << stage << " to " << val << std::endl;
	if (stage != 0) throw std::runtime_error("In a multiopt_scan, p2 should only be set at stage 0.");
	
      } else if (params[j].compare("B2c") == 0) {
	mo.opts[stage].q.B2c = val;
	if (verbose > 1) std::cout << "Setting B2c for opt stage " << stage << " to " << val << std::endl;
	if (stage != 0) throw std::runtime_error("In a multiopt_scan, B2c should only be set at stage 0.");

      } else if (params[j].compare("B2s") == 0) {
	mo.opts[stage].q.B2s = val;
	if (verbose > 1) std::cout << "Setting B2s for opt stage " << stage << " to " << val << std::endl;
	if (stage != 0) throw std::runtime_error("In a multiopt_scan, B2s should only be set at stage 0.");

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
  passed_filters = true;
  if (!keep_all) {
      
    if (mo.opts[index].q.grid_min_R0 < min_R0_to_keep) {
      filters_local[REJECTED_DUE_TO_R0]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to min_R0." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed min_R0 filter." << std::endl;
    
    if (std::abs(mo.opts[index].q.iota) < min_iota_to_keep) {
      filters_local[REJECTED_DUE_TO_IOTA]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to iota." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed iota filter." << std::endl;
    
    if (mo.opts[index].q.grid_max_elongation > max_elongation_to_keep) {
      filters_local[REJECTED_DUE_TO_ELONGATION]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to elongation. "
				 << mo.opts[index].q.grid_max_elongation << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed elongation filter. "
				      << mo.opts[index].q.grid_max_elongation << std::endl;
    
    if (mo.opts[index].q.grid_min_L_grad_B < min_L_grad_B_to_keep) {
      filters_local[REJECTED_DUE_TO_L_GRAD_B]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to L_grad_B." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed L_grad_B filter." << std::endl;
    
    if (mo.opts[index].q.B20_grid_variation > max_B20_variation_to_keep) {
      filters_local[REJECTED_DUE_TO_B20_VARIATION]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to B20 variation. "
				 << mo.opts[index].q.B20_grid_variation << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed B20 variation filter. "
				      << mo.opts[index].q.B20_grid_variation << std::endl;
    
    if (mo.opts[index].q.d2_volume_d_psi2 > max_d2_volume_d_psi2_to_keep) {
      filters_local[REJECTED_DUE_TO_D2_VOLUME_D_PSI2]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to d2_volume_d_psi2." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed d2_volume_d_psi2 filter." << std::endl;
    
    if (mo.opts[index].q.DMerc_times_r2 < min_DMerc_times_r2_to_keep) {
      filters_local[REJECTED_DUE_TO_DMERC]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to DMerc." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed DMerc filter." << std::endl;
    
    if (mo.opts[index].q.grid_min_L_grad_grad_B < min_L_grad_grad_B_to_keep) {
      filters_local[REJECTED_DUE_TO_L_GRAD_GRAD_B]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to L_grad_grad_B." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed L_grad_grad_B filter." << std::endl;
    
    if (mo.opts[index].q.r_singularity_robust < min_r_singularity_to_keep) {
      filters_local[REJECTED_DUE_TO_R_SINGULARITY]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to r_singularity. "
				 << mo.opts[index].q.r_singularity_robust << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed r_singularity filter. "
				      << mo.opts[index].q.r_singularity_robust << std::endl;
    
    if (mo.opts[index].q.grid_max_XY2 > max_XY2_to_keep) {
      filters_local[REJECTED_DUE_TO_MAX_XY2]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to max XY2." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed max XY2 filter." << std::endl;
    
    if (mo.opts[index].q.grid_max_Z2 > max_Z2_to_keep) {
      filters_local[REJECTED_DUE_TO_MAX_Z2]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to max Z2." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed max Z2 filter." << std::endl;
    
    if (mo.opts[index].q.grid_max_XY3 > max_XY3_to_keep) {
      filters_local[REJECTED_DUE_TO_MAX_XY3]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to max XY3." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed max XY3 filter." << std::endl;
    
    if (mo.opts[index].q.grid_max_d_XY2_d_varphi > max_d_XY2_d_varphi_to_keep) {
      filters_local[REJECTED_DUE_TO_MAX_D_XY2_D_VARPHI]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to max d XY2 / d varphi." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed max d XY2 / d varphi filter." << std::endl;
    
    if (mo.opts[index].q.grid_max_d_Z2_d_varphi > max_d_Z2_d_varphi_to_keep) {
      filters_local[REJECTED_DUE_TO_MAX_D_Z2_D_VARPHI]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to max d Z2 / d varphi." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed max d Z2 / d varphi filter." << std::endl;
    
    if (mo.opts[index].q.grid_max_d_XY3_d_varphi > max_d_XY3_d_varphi_to_keep) {
      filters_local[REJECTED_DUE_TO_MAX_D_XY3_D_VARPHI]++;
      if (verbose > 1) std::cout << "Rejecting this configuration due to max d XY3 / d varphi." << std::endl;
      passed_filters = false;	
    } else if (verbose > 1) std::cout << "Passed max d XY3 / d varphi filter." << std::endl;
    
  }
  if (passed_filters) filters_local[KEPT]++;

  // Package the results up for sending via MPI.
  // The order here must match the order in filter_global_arrays().
  int_parameters_single[0] = j_scan;
  int_parameters_single[1] = (int)passed_filters;
  int_parameters_single[2] = mo.n_evals;
  int_parameters_single[3] = mo.opts[index].q.helicity;
  for (j = 0; j < N_FILTERS; j++)
    int_parameters_single[j + n_int_parameters_base] = filters_local[j];
    
  end_time_single = std::chrono::steady_clock::now();
  elapsed = end_time_single - start_time_single;
  parameters_single[0] = elapsed.count();
  
  parameters_single[1 ] = mo.opts[index].q.eta_bar;
  parameters_single[2 ] = mo.opts[index].q.sigma0;
  parameters_single[3 ] = mo.opts[index].q.B2c;
  parameters_single[4 ] = mo.opts[index].q.B2s;
  parameters_single[5 ] = mo.opts[index].q.grid_min_R0;
  parameters_single[6 ] = mo.opts[index].q.grid_max_curvature;
  parameters_single[7 ] = mo.opts[index].q.iota;
  parameters_single[8 ] = mo.opts[index].q.grid_max_elongation;
  parameters_single[9 ] = mo.opts[index].q.grid_min_L_grad_B;
  parameters_single[10] = mo.opts[index].q.grid_min_L_grad_grad_B;
  parameters_single[11] = mo.opts[index].q.r_singularity_robust;
  parameters_single[12] = mo.opts[index].q.d2_volume_d_psi2;
  parameters_single[13] = mo.opts[index].q.DMerc_times_r2;
  parameters_single[14] = mo.opts[index].q.B20_grid_variation;
  parameters_single[15] = mo.opts[index].q.B20_residual;
  parameters_single[16] = mo.opts[index].q.B20_mean;
  parameters_single[17] = mo.opts[index].q.standard_deviation_of_R;
  parameters_single[18] = mo.opts[index].q.standard_deviation_of_Z;
  parameters_single[19] = mo.opts[index].q.grid_max_XY2;
  parameters_single[20] = mo.opts[index].q.grid_max_Z2;
  parameters_single[21] = mo.opts[index].q.grid_max_XY3;
  parameters_single[22] = mo.opts[index].q.grid_max_d_XY2_d_varphi;
  parameters_single[23] = mo.opts[index].q.grid_max_d_Z2_d_varphi;
  parameters_single[24] = mo.opts[index].q.grid_max_d_XY3_d_varphi;
  parameters_single[25] = mo.opts[index].q.grid_max_d2_XY2_d_varphi2;
  parameters_single[26] = mo.opts[index].q.grid_max_d2_XY3_d_varphi2;
  parameters_single[27] = mo.opts[index].q.axis_length;
  parameters_single[28] = mo.opts[index].q.p2;
  parameters_single[29] = mo.opts[0].iter_eta_bar[0];
  parameters_single[30] = mo.opts[0].iter_B2c[0];
  parameters_single[31] = mo.opts[0].iter_B2s[0];
    
  parameters_single[32] = mo.opts[index].weight_B20;
  parameters_single[33] = mo.opts[index].weight_iota;
  parameters_single[34] = mo.opts[index].target_iota;
  parameters_single[35] = mo.opts[index].weight_elongation;
  parameters_single[36] = mo.opts[index].weight_curvature;
  parameters_single[37] = mo.opts[index].weight_R0;
  parameters_single[38] = mo.opts[index].min_R0;
  parameters_single[39] = mo.opts[index].weight_d2_volume_d_psi2;
  parameters_single[40] = mo.opts[index].max_d2_volume_d_psi2;
  parameters_single[41] = mo.opts[index].weight_DMerc_times_r2;
  parameters_single[42] = mo.opts[index].min_DMerc_times_r2;
  parameters_single[43] = mo.opts[index].weight_XY2;
  parameters_single[44] = mo.opts[index].weight_XY2Prime;
  parameters_single[45] = mo.opts[index].weight_XY2PrimePrime;
  parameters_single[46] = mo.opts[index].weight_Z2;
  parameters_single[47] = mo.opts[index].weight_Z2Prime;
  parameters_single[48] = mo.opts[index].weight_XY3;
  parameters_single[49] = mo.opts[index].weight_XY3Prime;
  parameters_single[50] = mo.opts[index].weight_XY3PrimePrime;
  parameters_single[51] = mo.opts[index].weight_grad_B;
  parameters_single[52] = mo.opts[index].weight_grad_grad_B;
  parameters_single[53] = mo.opts[index].weight_r_singularity;
  parameters_single[54] = mo.opts[index].weight_axis_length;
  parameters_single[55] = mo.opts[index].target_axis_length;
  parameters_single[56] = mo.opts[index].weight_standard_deviation_of_R;
  parameters_single[57] = mo.opts[index].weight_B20_mean;

  for (j = 0; j < axis_nmax_plus_1; j++) {
    parameters_single[j + 0 * axis_nmax_plus_1 + n_parameters_base] = mo.opts[index].q.R0c[j];
    parameters_single[j + 1 * axis_nmax_plus_1 + n_parameters_base] = mo.opts[index].q.R0s[j];
    parameters_single[j + 2 * axis_nmax_plus_1 + n_parameters_base] = mo.opts[index].q.Z0c[j];
    parameters_single[j + 3 * axis_nmax_plus_1 + n_parameters_base] = mo.opts[index].q.Z0s[j];
  }
    
  end_time_single = std::chrono::steady_clock::now();
  elapsed = end_time_single - start_time_single;
  if (verbose > 0)
    std::cout << "Proc " << mpi_rank << " eval of index " << j_scan << " took " << elapsed.count() << " seconds" << std::endl;
}

void MultiOptScan::print_status() {
  int j;
  std::chrono::time_point<std::chrono::steady_clock> end_time;
  std::chrono::duration<double> elapsed;
  
  n_scan = filters[KEPT];
  
  end_time = std::chrono::steady_clock::now();
  elapsed = end_time - start_time;
  std::cout  << std::endl << "Status at " << elapsed.count() << " seconds:" << std::endl;
  qscfloat fraction_completed = ((qscfloat)filters[ATTEMPTS]) / n_scan_all;
  std::cout << filters[ATTEMPTS] << " of " << n_scan_all << " (" << fraction_completed << ") configs completed." << std::endl;
  std::cout << "Expected total scan time: " << elapsed.count() / fraction_completed << std::endl;
  std::cout << "CPU seconds for solves: " << total_cpu_seconds
	    << "  Wallclock time for that should be " << total_cpu_seconds / (n_procs - 1) << " seconds" << std::endl;
  std::cout << "Avg time per multiopt: " << total_cpu_seconds / filters[ATTEMPTS] << " seconds" << std::endl;
  std::cout << "Total number of function evaluations: " << n_evals
	    << "  Avg time per eval: " << total_cpu_seconds / n_evals << " seconds" << std::endl;
    
  std::cout << "Attempts on each proc:";
  for (j = 0; j < n_procs; j++) std::cout << " " << attempts_per_proc[j];
  std::cout << std::endl;
  std::cout << "# solves kept on each proc:";
  for (j = 0; j < n_procs; j++) std::cout << " " << n_solves_kept[j];
  std::cout << std::endl;

  for (j = 0; j < N_FILTERS; j++) filter_fractions[j] = ((qscfloat)filters[j]) / filters[0];

  int width = 13;
  std::cout << std::setprecision(4) << std::endl;
  std::cout << "Summary of scan results:                           (fractions in parentheses)" << std::endl;
  std::cout << "  Configurations attempted:          " << std::setw(width) << filters[ATTEMPTS] << std::endl;
  std::cout << "  Kept:                              " << std::setw(width) << n_scan
	    << " (" << filter_fractions[KEPT] << ")" << std::endl;
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
  std::cout << "  Rejected due to max XY2:           " << std::setw(width) << filters[REJECTED_DUE_TO_MAX_XY2]
	    << " (" << filter_fractions[REJECTED_DUE_TO_MAX_XY2] << ")" << std::endl;
  std::cout << "  Rejected due to max Z2:            " << std::setw(width) << filters[REJECTED_DUE_TO_MAX_Z2]
	    << " (" << filter_fractions[REJECTED_DUE_TO_MAX_Z2] << ")" << std::endl;
  std::cout << "  Rejected due to max XY3:           " << std::setw(width) << filters[REJECTED_DUE_TO_MAX_XY3]
	    << " (" << filter_fractions[REJECTED_DUE_TO_MAX_XY3] << ")" << std::endl;
  std::cout << "  Rejected due to max dXY2/dvarphi:  " << std::setw(width) << filters[REJECTED_DUE_TO_MAX_D_XY2_D_VARPHI]
	    << " (" << filter_fractions[REJECTED_DUE_TO_MAX_D_XY2_D_VARPHI] << ")" << std::endl;
  std::cout << "  Rejected due to max dZ2/dvarphi:   " << std::setw(width) << filters[REJECTED_DUE_TO_MAX_D_Z2_D_VARPHI]
	    << " (" << filter_fractions[REJECTED_DUE_TO_MAX_D_Z2_D_VARPHI] << ")" << std::endl;
  std::cout << "  Rejected due to max dXY3/dvarphi:  " << std::setw(width) << filters[REJECTED_DUE_TO_MAX_D_XY3_D_VARPHI]
	    << " (" << filter_fractions[REJECTED_DUE_TO_MAX_D_XY3_D_VARPHI] << ")" << std::endl;
  std::cout << std::endl;
}

void MultiOptScan::filter_global_arrays() {
  int j, k, j_global;
  std::chrono::time_point<std::chrono::steady_clock> start_time_filter, end_time;
  std::chrono::duration<double> elapsed;
  
  start_time_filter = std::chrono::steady_clock::now();
  
  n_scan = filters[KEPT];
  
  scan_helicity.resize(n_scan, 0);
  scan_n_evals.resize(n_scan, 0);

  scan_eta_bar.resize(n_scan, 0.0);
  scan_sigma0.resize(n_scan, 0.0);
  scan_B2c.resize(n_scan, 0.0);
  scan_B2s.resize(n_scan, 0.0);
  scan_min_R0.resize(n_scan, 0.0);
  scan_max_curvature.resize(n_scan, 0.0);
  scan_iota.resize(n_scan, 0.0);
  scan_max_elongation.resize(n_scan, 0.0);
  scan_min_L_grad_B.resize(n_scan, 0.0);
  scan_min_L_grad_grad_B.resize(n_scan, 0.0);
  scan_r_singularity.resize(n_scan, 0.0);
  scan_d2_volume_d_psi2.resize(n_scan, 0.0);
  scan_DMerc_times_r2.resize(n_scan, 0.0);
  scan_B20_variation.resize(n_scan, 0.0);
  scan_B20_residual.resize(n_scan, 0.0);
  scan_B20_mean.resize(n_scan, 0.0);
  scan_standard_deviation_of_R.resize(n_scan, 0.0);
  scan_standard_deviation_of_Z.resize(n_scan, 0.0);
  
  scan_max_XY2.resize(n_scan, 0.0);
  scan_max_Z2.resize(n_scan, 0.0);
  scan_max_XY3.resize(n_scan, 0.0);
  scan_max_d_XY2_d_varphi.resize(n_scan, 0.0);
  scan_max_d_Z2_d_varphi.resize(n_scan, 0.0);
  scan_max_d_XY3_d_varphi.resize(n_scan, 0.0);
  scan_max_d2_XY2_d_varphi2.resize(n_scan, 0.0);
  scan_max_d2_XY3_d_varphi2.resize(n_scan, 0.0);
  scan_axis_length.resize(n_scan, 0.0);
  scan_p2.resize(n_scan, 0.0);
  scan_initial_eta_bar.resize(n_scan, 0.0);
  scan_initial_B2c.resize(n_scan, 0.0);
  scan_initial_B2s.resize(n_scan, 0.0);

  scan_R0c.resize(axis_nmax_plus_1, n_scan, 0.0);
  scan_R0s.resize(axis_nmax_plus_1, n_scan, 0.0);
  scan_Z0c.resize(axis_nmax_plus_1, n_scan, 0.0);
  scan_Z0s.resize(axis_nmax_plus_1, n_scan, 0.0);

  scan_weight_B20.resize(n_scan, 0.0);
  scan_weight_iota.resize(n_scan, 0.0);
  scan_target_iota.resize(n_scan, 0.0);
  scan_weight_elongation.resize(n_scan, 0.0);
  scan_weight_curvature.resize(n_scan, 0.0);
  scan_weight_R0.resize(n_scan, 0.0);
  scan_target_min_R0.resize(n_scan, 0.0);
  scan_weight_d2_volume_d_psi2.resize(n_scan, 0.0);
  scan_max_d2_volume_d_psi2.resize(n_scan, 0.0);
  scan_weight_DMerc_times_r2.resize(n_scan, 0.0);
  scan_min_DMerc_times_r2.resize(n_scan, 0.0);
  scan_weight_XY2.resize(n_scan, 0.0);
  scan_weight_XY2Prime.resize(n_scan, 0.0);
  scan_weight_XY2PrimePrime.resize(n_scan, 0.0);
  scan_weight_Z2.resize(n_scan, 0.0);
  scan_weight_Z2Prime.resize(n_scan, 0.0);
  scan_weight_XY3.resize(n_scan, 0.0);
  scan_weight_XY3Prime.resize(n_scan, 0.0);
  scan_weight_XY3PrimePrime.resize(n_scan, 0.0);
  scan_weight_grad_B.resize(n_scan, 0.0);
  scan_weight_grad_grad_B.resize(n_scan, 0.0);
  scan_weight_r_singularity.resize(n_scan, 0.0);
  scan_weight_axis_length.resize(n_scan, 0.0);
  scan_target_axis_length.resize(n_scan, 0.0);
  scan_weight_standard_deviation_of_R.resize(n_scan, 0.0);
  scan_weight_B20_mean.resize(n_scan, 0.0);
    
  // Unpack parameters.
  // The order here must match the order at the end of eval_scan_index()!
  j = -1;
  for (j_global = 0; j_global < n_scan_all; j_global++) {
    if (int_parameters[1 + j_global * n_int_parameters_base] < 1) continue;
    j++;
    
    scan_n_evals[j]  = int_parameters[2 + j_global * n_int_parameters_base];
    scan_helicity[j] = int_parameters[3 + j_global * n_int_parameters_base];
      
    // run time is #0
    scan_eta_bar[j]           = parameters( 1, j_global);
    scan_sigma0[j]            = parameters( 2, j_global);
    scan_B2c[j]               = parameters( 3, j_global);
    scan_B2s[j]               = parameters( 4, j_global);
    scan_min_R0[j]            = parameters( 5, j_global);
    scan_max_curvature[j]     = parameters( 6, j_global);
    scan_iota[j]              = parameters( 7, j_global);
    scan_max_elongation[j]    = parameters( 8, j_global);
    scan_min_L_grad_B[j]      = parameters( 9, j_global);
    scan_min_L_grad_grad_B[j] = parameters(10, j_global);
    scan_r_singularity[j]     = parameters(11, j_global);
    scan_d2_volume_d_psi2[j]  = parameters(12, j_global);
    scan_DMerc_times_r2[j]    = parameters(13, j_global);
    scan_B20_variation[j]     = parameters(14, j_global);
    scan_B20_residual[j]      = parameters(15, j_global);
    scan_B20_mean[j]          = parameters(16, j_global);
    scan_standard_deviation_of_R[j] = parameters(17, j_global);
    scan_standard_deviation_of_Z[j] = parameters(18, j_global);
    scan_max_XY2[j]            = parameters(19, j_global);
    scan_max_Z2[j]             = parameters(20, j_global);
    scan_max_XY3[j]            = parameters(21, j_global);
    scan_max_d_XY2_d_varphi[j] = parameters(22, j_global);
    scan_max_d_Z2_d_varphi[j]  = parameters(23, j_global);
    scan_max_d_XY3_d_varphi[j] = parameters(24, j_global);
    scan_max_d2_XY2_d_varphi2[j] = parameters(25, j_global);
    scan_max_d2_XY3_d_varphi2[j] = parameters(26, j_global);
    scan_axis_length[j]        = parameters(27, j_global);
    scan_p2[j]                 = parameters(28, j_global);
    scan_initial_eta_bar[j]    = parameters(29, j_global);
    scan_initial_B2c[j]        = parameters(30, j_global);
    scan_initial_B2s[j]        = parameters(31, j_global);
    
    scan_weight_B20[j]                     = parameters(32, j_global);
    scan_weight_iota[j]                    = parameters(33, j_global);
    scan_target_iota[j]                    = parameters(34, j_global);
    scan_weight_elongation[j]              = parameters(35, j_global);
    scan_weight_curvature[j]               = parameters(36, j_global);
    scan_weight_R0[j]                      = parameters(37, j_global);
    scan_target_min_R0[j]                  = parameters(38, j_global);
    scan_weight_d2_volume_d_psi2[j]        = parameters(39, j_global);
    scan_max_d2_volume_d_psi2[j]           = parameters(40, j_global);
    scan_weight_DMerc_times_r2[j]          = parameters(41, j_global);
    scan_min_DMerc_times_r2[j]             = parameters(42, j_global);
    scan_weight_XY2[j]                     = parameters(43, j_global);
    scan_weight_XY2Prime[j]                = parameters(44, j_global);
    scan_weight_XY2PrimePrime[j]           = parameters(45, j_global);
    scan_weight_Z2[j]                      = parameters(46, j_global);
    scan_weight_Z2Prime[j]                 = parameters(47, j_global);
    scan_weight_XY3[j]                     = parameters(48, j_global);
    scan_weight_XY3Prime[j]                = parameters(49, j_global);
    scan_weight_XY3PrimePrime[j]           = parameters(50, j_global);
    scan_weight_grad_B[j]                  = parameters(51, j_global);
    scan_weight_grad_grad_B[j]             = parameters(52, j_global);
    scan_weight_r_singularity[j]           = parameters(53, j_global);
    scan_weight_axis_length[j]             = parameters(54, j_global);
    scan_target_axis_length[j]             = parameters(55, j_global);
    scan_weight_standard_deviation_of_R[j] = parameters(56, j_global);
    scan_weight_B20_mean[j]                = parameters(57, j_global);
    
    for (k = 0; k < axis_nmax_plus_1; k++) {
      scan_R0c(k, j) = parameters(k + 0 * axis_nmax_plus_1 + n_parameters_base, j_global);
      scan_R0s(k, j) = parameters(k + 1 * axis_nmax_plus_1 + n_parameters_base, j_global);
      scan_Z0c(k, j) = parameters(k + 2 * axis_nmax_plus_1 + n_parameters_base, j_global);
      scan_Z0s(k, j) = parameters(k + 3 * axis_nmax_plus_1 + n_parameters_base, j_global);
    }
  }
  if (j + 1 != n_scan) {
    std::cout << "Error! mismatch in number of configs saved. n_scan=" << n_scan << " j=" << j << std::endl;
    throw std::runtime_error("mismatch in number of configs saved.");
  }
  end_time = std::chrono::steady_clock::now();
  elapsed = end_time - start_time_filter;
  std::cout << "Time for filter_global_arrays: " << elapsed.count() << " seconds:" << std::endl;  
}
