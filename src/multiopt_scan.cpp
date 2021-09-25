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
  save_period = 10;
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
  min_DMerc_to_keep = -1.0e+30;
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
  MPI_Status mpi_status;
  char str_buffer[100];
  
  // Initialize MPI-related variables
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Comm_size(mpi_comm, &n_procs);
  if (n_procs < 2)
    throw std::runtime_error("For MultiOptScan, the number of MPI processes must be at least 2.");
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
    if (verbose > 0 && proc0)
      std::cout << "Values for parameter " << params[j] << ": " << params_vals[j] << std::endl;
  }
  if (verbose > 0 && proc0) std::cout << "Total number of points in scan: " << n_scan_all << std::endl;

  for (j = 0; j < N_FILTERS; j++) filters_local[j] = 0;
  axis_nmax_plus_1 = mo_ref.opts[0].q.R0c.size();
  for (j = 0; j < mo_ref.opts.size(); j++) axis_nmax_plus_1 += mo_ref.opts[j].fourier_refine;

  if (n_procs - 1 > n_scan_all)
    throw std::runtime_error("For MultiOptScan, the number of MPI processes cannot exceed n_scan_all + 1.");
  
  if (n_scan_all >= n_procs) {
    scan_index_min = (mpi_rank * n_scan_all) / n_procs;
    scan_index_max = ((mpi_rank + 1) * n_scan_all) / n_procs - 1;
  } else {
    // An uncommon case: There are more procs than solves to do                                    
    scan_index_min = std::min((big)mpi_rank, n_scan_all);
    if (mpi_rank < n_scan_all) {
      scan_index_max = scan_index_min;
    } else {
      scan_index_max = scan_index_min - 1;
    }
  }
  n_scan_local = scan_index_max - scan_index_min + 1;

  // Print the processor assignments in a coordinated manner.
  std::string proc_assignments_string;
  std::stringstream ss;
  ss << "Proc" << std::setw(5) << mpi_rank << " of" << std::setw(5) << n_procs << " will handle points"
     << std::setw(9) << scan_index_min << " to" << std::setw(9) << scan_index_max
     << " (n =" << std::setw(9) << n_scan_local << ")";
  proc_assignments_string = ss.str();
  //std::cout << " Length: " << proc_assignments_string.length() << std::endl;
  int str_length = 72;
  assert (proc_assignments_string.length() == str_length);
  if (proc0) {
    std::cout << proc_assignments_string << std::endl;
    for (j = 1; j < n_procs; j++) {
      MPI_Recv(str_buffer, str_length, MPI_CHAR, j, j, mpi_comm, &mpi_status);
      std::cout << std::string(str_buffer) << std::endl;
    }
  } else {
    MPI_Send(proc_assignments_string.data(), str_length, MPI_CHAR, 0, mpi_rank, mpi_comm);
  }

  // Begin test code for MPI master-slave system.
  int intbuf[1];
  int proc_that_finished;
  if (proc0) {
    // Send off the first batch of work:
    for (j = 1; j < n_procs; j++) {
      intbuf[0] = j - 1;
      MPI_Send(intbuf, 1, MPI_INT, j, 0, mpi_comm);
    }
    // Handle the rest of the tasks:
    for (j = n_procs - 1; j < n_scan_all; j++) {
      // Wait for any proc to finish:
      MPI_Recv(intbuf, 1, MPI_INT, MPI_ANY_SOURCE, 0, mpi_comm, &mpi_status);
      proc_that_finished = intbuf[0];
      std::cout << "proc0: I see that proc " << proc_that_finished << " just finished, so I'll send him index " << j << std::endl;
      intbuf[0] = j;
      // Give that proc the next piece of work:
      MPI_Send(intbuf, 1, MPI_INT, proc_that_finished, 0, mpi_comm);
    }
    // Receive the final completion messages from all other procs:
    for (j = 1; j < n_procs; j++) {
      MPI_Recv(intbuf, 1, MPI_INT, MPI_ANY_SOURCE, 0, mpi_comm, &mpi_status);
      proc_that_finished = intbuf[0];
      std::cout << "proc0: I see that proc " << proc_that_finished << " just finished. No further work." << std::endl;
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
	std::cout << "proc " << mpi_rank << " is exiting." << std::endl;
      } else {
	std::cout << "proc " << mpi_rank << " is handling index " << intbuf[0] << std::endl;
	// Notify proc0 that we are done.
	intbuf[0] = mpi_rank;
	MPI_Send(intbuf, 1, MPI_INT, 0, 0, mpi_comm);
      }
    }
  }
}

void MultiOptScan::scan() {
  const int n_parameters = 37;
  const int n_int_parameters = 1;
  const int n_fourier_parameters = axis_nmax_plus_1 * 4;
  Matrix parameters_local(n_parameters, n_scan_local);
  Matrix fourier_parameters_local(n_fourier_parameters, n_scan_local);
  std::valarray<int> int_parameters_local(n_int_parameters * n_scan_local);
  big j_scan, j_scan_global, jmod;
  int j, k, stage, stage_min, stage_max, index;
  std::valarray<int> indices;
  qscfloat val;
  std::chrono::time_point<std::chrono::steady_clock> end_time;
  std::chrono::duration<double> elapsed;
  
  // Some indices of interest:
  // j_scan_global is the index into the complete tensor-product set of points.
  // j_scan is the index into the current saved solution on this process, after filtering out some configs.
  // filters_local[ATTEMPTS] is the number of configurations considered so far on this process.
  
  indices.resize(ndim);
  j_scan = 0;
  start_time = std::chrono::steady_clock::now();
  for (j_scan_global = scan_index_min; j_scan_global <= scan_index_max; j_scan_global++) {
    filters_local[ATTEMPTS]++;

    // Don't collect results if we are close to finishing (reaching
    // n_scan_all / n_procs attempts), since different procs may have
    // n_scan_local values differing by 1. Be careful to put +2 on the
    // LHS instead of -2 on the RHS since bigs are unsigned.
    if ((filters_local[ATTEMPTS] % save_period == 0) && (filters_local[ATTEMPTS] + 2 < n_scan_all / n_procs)) {
      end_time = std::chrono::steady_clock::now();
      elapsed = end_time - start_time;
      std::cout << "proc " << mpi_rank << " calling collect_results in loop after "
		<< elapsed.count() << " seconds" << std::endl;
      collect_results(n_parameters, parameters_local, fourier_parameters_local,
                      n_int_parameters, int_parameters_local, j_scan);
    }
    
    // From the global scan index j_scan_global, determine the indices for each dimension of the scan:
    jmod = j_scan_global;
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
	filters_local[REJECTED_DUE_TO_R0]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to min_R0." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed min_R0 filter." << std::endl;
      
      if (std::abs(mo.opts[index].q.iota) < min_iota_to_keep) {
	filters_local[REJECTED_DUE_TO_IOTA]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to iota." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed iota filter." << std::endl;
      
      if (mo.opts[index].q.grid_max_elongation > max_elongation_to_keep) {
	filters_local[REJECTED_DUE_TO_ELONGATION]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to elongation." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed elongation filter." << std::endl;
      
      if (mo.opts[index].q.grid_min_L_grad_B < min_L_grad_B_to_keep) {
	filters_local[REJECTED_DUE_TO_L_GRAD_B]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to L_grad_B." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed L_grad_B filter." << std::endl;
      
      if (mo.opts[index].q.B20_grid_variation > max_B20_variation_to_keep) {
	filters_local[REJECTED_DUE_TO_B20_VARIATION]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to B20 variation." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed B20 variation filter." << std::endl;
      
      if (mo.opts[index].q.d2_volume_d_psi2 > max_d2_volume_d_psi2_to_keep) {
	filters_local[REJECTED_DUE_TO_D2_VOLUME_D_PSI2]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to d2_volume_d_psi2." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed d2_volume_d_psi2 filter." << std::endl;
      
      if (mo.opts[index].q.DMerc_times_r2 < min_DMerc_to_keep) {
	filters_local[REJECTED_DUE_TO_DMERC]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to DMerc." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed DMerc filter." << std::endl;
      
      if (mo.opts[index].q.grid_min_L_grad_grad_B < min_L_grad_grad_B_to_keep) {
	filters_local[REJECTED_DUE_TO_L_GRAD_GRAD_B]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to L_grad_grad_B." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed L_grad_grad_B filter." << std::endl;
      
      if (mo.opts[index].q.r_singularity_robust < min_r_singularity_to_keep) {
	filters_local[REJECTED_DUE_TO_R_SINGULARITY]++;
	if (verbose > 1) std::cout << "Rejecting this configuration due to r_singularity." << std::endl;
	continue;	
      } else if (verbose > 1) std::cout << "Passed r_singularity filter." << std::endl;
      
    }

    // If we make it this far, then save the configuration.
    parameters_local(0 , j_scan) = mo.opts[index].q.eta_bar;
    parameters_local(1 , j_scan) = mo.opts[index].q.sigma0;
    parameters_local(2 , j_scan) = mo.opts[index].q.B2c;
    parameters_local(3 , j_scan) = mo.opts[index].q.B2s;
    parameters_local(4 , j_scan) = mo.opts[index].q.grid_min_R0;
    parameters_local(5 , j_scan) = mo.opts[index].q.grid_max_curvature;
    parameters_local(6 , j_scan) = mo.opts[index].q.iota;
    parameters_local(7 , j_scan) = mo.opts[index].q.grid_max_elongation;
    parameters_local(8 , j_scan) = mo.opts[index].q.grid_min_L_grad_B;
    parameters_local(9 , j_scan) = mo.opts[index].q.grid_min_L_grad_grad_B;
    parameters_local(10, j_scan) = mo.opts[index].q.r_singularity_robust;
    parameters_local(11, j_scan) = mo.opts[index].q.d2_volume_d_psi2;
    parameters_local(12, j_scan) = mo.opts[index].q.DMerc_times_r2;
    parameters_local(13, j_scan) = mo.opts[index].q.B20_grid_variation;
    parameters_local(14, j_scan) = mo.opts[index].q.B20_residual;
    parameters_local(15, j_scan) = mo.opts[index].q.standard_deviation_of_R;
    parameters_local(16, j_scan) = mo.opts[index].q.standard_deviation_of_Z;
    
    parameters_local(17, j_scan) = mo.opts[index].weight_B20;
    parameters_local(18, j_scan) = mo.opts[index].weight_iota;
    parameters_local(19, j_scan) = mo.opts[index].target_iota;
    parameters_local(20, j_scan) = mo.opts[index].weight_elongation;
    parameters_local(21, j_scan) = mo.opts[index].weight_curvature;
    parameters_local(22, j_scan) = mo.opts[index].weight_R0;
    parameters_local(23, j_scan) = mo.opts[index].min_R0;
    parameters_local(24, j_scan) = mo.opts[index].weight_d2_volume_d_psi2;
    parameters_local(25, j_scan) = mo.opts[index].max_d2_volume_d_psi2;
    parameters_local(26, j_scan) = mo.opts[index].weight_XY2;
    parameters_local(27, j_scan) = mo.opts[index].weight_XY2Prime;
    parameters_local(28, j_scan) = mo.opts[index].weight_Z2;
    parameters_local(29, j_scan) = mo.opts[index].weight_Z2Prime;
    parameters_local(30, j_scan) = mo.opts[index].weight_XY3;
    parameters_local(31, j_scan) = mo.opts[index].weight_XY3Prime;
    parameters_local(32, j_scan) = mo.opts[index].weight_grad_B;
    parameters_local(33, j_scan) = mo.opts[index].weight_grad_grad_B;
    parameters_local(34, j_scan) = mo.opts[index].weight_r_singularity;
    parameters_local(35, j_scan) = mo.opts[index].weight_axis_length;
    parameters_local(36, j_scan) = mo.opts[index].weight_standard_deviation_of_R;

    int_parameters_local[0 + n_int_parameters * j_scan] = mo.opts[index].q.helicity;
    
    for (j = 0; j < axis_nmax_plus_1; j++) {
      fourier_parameters_local(j + 0 * axis_nmax_plus_1, j_scan) = mo.opts[index].q.R0c[j];
      fourier_parameters_local(j + 1 * axis_nmax_plus_1, j_scan) = mo.opts[index].q.R0s[j];
      fourier_parameters_local(j + 2 * axis_nmax_plus_1, j_scan) = mo.opts[index].q.Z0c[j];
      fourier_parameters_local(j + 3 * axis_nmax_plus_1, j_scan) = mo.opts[index].q.Z0s[j];
    }
    
    j_scan++;

  }

  end_time = std::chrono::steady_clock::now();
  elapsed = end_time - start_time;
  std::cout << "Proc " << mpi_rank << " finished after " << elapsed.count() << " seconds" << std::endl;

  // Collect final results:
  std::cout << "proc " << mpi_rank << " calling collect_results at end" << std::endl;
  collect_results(n_parameters, parameters_local, fourier_parameters_local,
                  n_int_parameters, int_parameters_local, j_scan);
}
