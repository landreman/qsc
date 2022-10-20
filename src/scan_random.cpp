#include <chrono>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "qsc.hpp"
#include "scan.hpp"
#include "random.hpp"

using namespace qsc;

void Scan::random() {
  const int n_parameters = 19;
  const int n_int_parameters = 1;
  const int axis_nmax_plus_1 = R0c_max.size();
  const int n_fourier_parameters = axis_nmax_plus_1 * 4;
  Matrix parameters_local(n_parameters, max_keep_per_proc);
  Matrix fourier_parameters_local(n_fourier_parameters, max_keep_per_proc);
  std::valarray<int> int_parameters_local(n_int_parameters * max_keep_per_proc);
  big j_scan = 0;
  int j, k;
  qscfloat R0_at_0, R0_at_half_period, val;
  int mpi_rank, n_procs;
  MPI_Status mpi_status;
  
  // Initialize MPI
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Comm_size(mpi_comm, &n_procs);
  bool proc0 = (mpi_rank == 0);
  
  // Initialize random distributions:
  Random random_eta_bar(deterministic, eta_bar_scan_option, eta_bar_min, eta_bar_max);
  Random random_sigma0(deterministic, sigma0_scan_option, sigma0_min, sigma0_max);
  Random random_B2c(deterministic, B2c_scan_option, B2c_min, B2c_max);
  Random random_B2s(deterministic, B2s_scan_option, B2s_min, B2s_max);
  Random random_p2(deterministic, p2_scan_option, p2_min, p2_max);
  random_eta_bar.set_to_nth(mpi_rank * max_attempts_per_proc + 0);
  random_sigma0.set_to_nth(mpi_rank * max_attempts_per_proc + 1);
  random_B2c.set_to_nth(mpi_rank * max_attempts_per_proc + 2);
  random_B2s.set_to_nth(mpi_rank * max_attempts_per_proc + 3);
  random_p2.set_to_nth(mpi_rank * max_attempts_per_proc + 4);
  Random* random_R0c[axis_nmax_plus_1];
  Random* random_R0s[axis_nmax_plus_1];
  Random* random_Z0c[axis_nmax_plus_1];
  Random* random_Z0s[axis_nmax_plus_1];
  for (j = 0; j < axis_nmax_plus_1; j++) {
    // For fourier_scan_option = "2 sided log", handle R0c[0] separately so it is not negative
    if (j == 0 && fourier_scan_option.compare(RANDOM_OPTION_2_SIDED_LOG) == 0) {
      random_R0c[j] = new Random(deterministic, RANDOM_OPTION_LOG, R0c_min[j], R0c_max[j]);
    } else {
      random_R0c[j] = new Random(deterministic, fourier_scan_option, R0c_min[j], R0c_max[j]);
    }
    random_R0s[j] = new Random(deterministic, fourier_scan_option, R0s_min[j], R0s_max[j]);
    random_Z0c[j] = new Random(deterministic, fourier_scan_option, Z0c_min[j], Z0c_max[j]);
    random_Z0s[j] = new Random(deterministic, fourier_scan_option, Z0s_min[j], Z0s_max[j]);
    random_R0c[j]->set_to_nth(mpi_rank * max_attempts_per_proc + 4);
    random_R0s[j]->set_to_nth(mpi_rank * max_attempts_per_proc + 5);
    random_Z0c[j]->set_to_nth(mpi_rank * max_attempts_per_proc + 6);
    random_Z0s[j]->set_to_nth(mpi_rank * max_attempts_per_proc + 7);
    // The set_to_nth() calls are so we can test that results are
    // independent of the number of MPI procs for deterministic
    // runs. Also, the +1, +2, ... in these calls is so the parameters
    // are a bit less correlated for deterministic runs.
  }

  for (j = 0; j < N_FILTERS; j++) filters_local[j] = 0;

  q.R0c.resize(axis_nmax_plus_1, 0.0);
  q.R0s.resize(axis_nmax_plus_1, 0.0);
  q.Z0c.resize(axis_nmax_plus_1, 0.0);
  q.Z0s.resize(axis_nmax_plus_1, 0.0);

  std::chrono::time_point<std::chrono::steady_clock> end_time, checkpoint_time;
  std::chrono::time_point<std::chrono::steady_clock> section_start_time, section_end_time;
  start_time = std::chrono::steady_clock::now();
  checkpoint_time = start_time;
  std::chrono::duration<double> elapsed;
  for (j = 0; j < N_TIMES; j++) timing_local[j] = 0.0;

  // Initialize the Qsc object:
  q.validate();
  q.allocate();
  
  bool keep_going = true;
  while (keep_going) {
    end_time = std::chrono::steady_clock::now();
    elapsed = end_time - start_time;
    if (elapsed.count() > max_seconds) keep_going = false;
    
    end_time = std::chrono::steady_clock::now();
    elapsed = end_time - checkpoint_time;
    if (elapsed.count() > save_period) {
      checkpoint_time = end_time;
      collect_results(n_parameters, parameters_local, fourier_parameters_local,
		      n_int_parameters, int_parameters_local, j_scan);
      write_netcdf();
    }
    
    if (max_attempts_per_proc > 0 && filters_local[ATTEMPTS] >= max_attempts_per_proc) break;
    
    filters_local[ATTEMPTS]++;

    // Pick random parameters. A small amount of time could be saved
    // if random numbers were requested later, only when needed, if
    // you make it past initial filters. However this makes it hard to
    // test that the results are independent of the # of MPI
    // processes, because then proc j would need a seed that depends
    // on how many cases pass the filters on proc j-1. I don't think
    // the random number generation is likely to take much of the
    // overall time for a scan, so let's just get the random values
    // here.
    section_start_time = std::chrono::steady_clock::now();
    q.eta_bar = random_eta_bar.get();
    q.sigma0 = random_sigma0.get();
    if (q.at_least_order_r2) {
      q.B2c = random_B2c.get();
      q.B2s = random_B2s.get();
      q.p2 = random_p2.get();
    }
    // Initialize axis, and do a crude check of whether R0 goes negative:
    R0_at_half_period = 0;
    for (j = 0; j < axis_nmax_plus_1; j++) {
      q.R0s[j] = random_R0s[j]->get();
      q.Z0s[j] = random_Z0s[j]->get();
      q.Z0c[j] = random_Z0c[j]->get();
      val = random_R0c[j]->get();
      q.R0c[j] = val;

      if (j % 2 == 0) {
	R0_at_half_period += val;
      } else {
	R0_at_half_period -= val;
      }
    }
    R0_at_0 = q.R0c.sum();
    section_end_time = std::chrono::steady_clock::now();
    elapsed = section_end_time - section_start_time;
    timing_local[TIME_RANDOM] += elapsed.count();
    if (R0_at_0 <= 0 || R0_at_half_period <= 0) {
      filters_local[REJECTED_DUE_TO_R0_CRUDE]++;
      continue;
    }

    section_start_time = std::chrono::steady_clock::now();
    q.init_axis();
    section_end_time = std::chrono::steady_clock::now();
    elapsed = section_end_time - section_start_time;
    timing_local[TIME_INIT_AXIS] += elapsed.count();
    if (!keep_all && q.grid_min_R0 < min_R0_to_keep) {
      filters_local[REJECTED_DUE_TO_R0]++;
      continue;
    }
    if (!keep_all) {
      // If helicity_to_keep == 0, only helicity == 0 will be kept.
      // If helicity_to_keep == 1, only helicity == 1 or -1 will be kept.
      // For any other value of helicity_to_keep, all helicity values will be kept.
      if ((helicity_to_keep == 0 && q.helicity != 0) || (helicity_to_keep == 1 && std::abs(q.helicity) != 1)) {
	filters_local[REJECTED_DUE_TO_HELICITY]++;
	continue;
      }
    }
    if (!keep_all && 1.0 / q.grid_max_curvature < min_L_grad_B_to_keep) {
      filters_local[REJECTED_DUE_TO_CURVATURE]++;
      continue;
    }

    // Here is the main O(r^1) solve:
    section_start_time = std::chrono::steady_clock::now();
    q.solve_sigma_equation();
    section_end_time = std::chrono::steady_clock::now();
    elapsed = section_end_time - section_start_time;
    timing_local[TIME_SIGMA_EQUATION] += elapsed.count();
    filters_local[N_SIGMA_EQ_SOLVES]++;
    
    section_start_time = std::chrono::steady_clock::now();
    q.r1_diagnostics();
    section_end_time = std::chrono::steady_clock::now();
    elapsed = section_end_time - section_start_time;
    timing_local[TIME_R1_DIAGNOSTICS] += elapsed.count();
    if (!keep_all && std::abs(q.iota) < min_iota_to_keep) {
      filters_local[REJECTED_DUE_TO_IOTA]++;
      continue;
    }
    if (!keep_all && q.grid_max_elongation > max_elongation_to_keep) {
      filters_local[REJECTED_DUE_TO_ELONGATION]++;
      continue;
    }
    if (!keep_all && q.grid_min_L_grad_B < min_L_grad_B_to_keep) {
      filters_local[REJECTED_DUE_TO_L_GRAD_B]++;
      continue;
    }

    if (q.at_least_order_r2) {
      // Here is the main O(r^2) solve:
      section_start_time = std::chrono::steady_clock::now();
      q.calculate_r2();
      filters_local[N_R2_SOLVES]++;
      section_end_time = std::chrono::steady_clock::now();
      elapsed = section_end_time - section_start_time;
      timing_local[TIME_CALCULATE_R2] += elapsed.count();

      // Filter results:
      if (!keep_all && q.B20_grid_variation > max_B20_variation_to_keep) {
	filters_local[REJECTED_DUE_TO_B20_VARIATION]++;
	continue;
      }

      section_start_time = std::chrono::steady_clock::now();
      q.mercier();
      section_end_time = std::chrono::steady_clock::now();
      elapsed = section_end_time - section_start_time;
      timing_local[TIME_MERCIER] += elapsed.count();
      if (!keep_all && q.d2_volume_d_psi2 > max_d2_volume_d_psi2_to_keep) {
	filters_local[REJECTED_DUE_TO_D2_VOLUME_D_PSI2]++;
	continue;
      }
      if (!keep_all && q.DMerc_times_r2 < min_DMerc_times_r2_to_keep) {
	filters_local[REJECTED_DUE_TO_DMERC]++;
	continue;
      }

      section_start_time = std::chrono::steady_clock::now();
      q.calculate_grad_grad_B_tensor();
      section_end_time = std::chrono::steady_clock::now();
      elapsed = section_end_time - section_start_time;
      timing_local[TIME_GRAD_GRAD_B_TENSOR] += elapsed.count();
      if (!keep_all && q.grid_min_L_grad_grad_B < min_L_grad_grad_B_to_keep) {
	filters_local[REJECTED_DUE_TO_L_GRAD_GRAD_B]++;
	continue;
      }

      section_start_time = std::chrono::steady_clock::now();
      q.calculate_r_singularity();
      section_end_time = std::chrono::steady_clock::now();
      elapsed = section_end_time - section_start_time;
      timing_local[TIME_R_SINGULARITY] += elapsed.count();
      if (!keep_all && q.r_singularity_robust < min_r_singularity_to_keep) {
	filters_local[REJECTED_DUE_TO_R_SINGULARITY]++;
	continue;
      }
      if (!keep_all && q.beta < min_beta_to_keep) {
	filters_local[REJECTED_DUE_TO_BETA]++;
	continue;
      }
    } // if at_least_order_r2

    // If we made it this far, then we found a keeper.
    parameters_local(0 , j_scan) = q.eta_bar;
    parameters_local(1 , j_scan) = q.sigma0;
    parameters_local(2 , j_scan) = q.B2c;
    parameters_local(3 , j_scan) = q.B2s;
    parameters_local(4 , j_scan) = q.p2;
    parameters_local(5 , j_scan) = q.grid_min_R0;
    parameters_local(6 , j_scan) = q.grid_max_curvature;
    parameters_local(7 , j_scan) = q.iota;
    parameters_local(8 , j_scan) = q.grid_max_elongation;
    parameters_local(9 , j_scan) = q.grid_min_L_grad_B;
    parameters_local(10, j_scan) = q.grid_min_L_grad_grad_B;
    parameters_local(11, j_scan) = q.r_singularity_robust;
    parameters_local(12, j_scan) = q.d2_volume_d_psi2;
    parameters_local(13, j_scan) = q.DMerc_times_r2;
    parameters_local(14, j_scan) = q.B20_grid_variation;
    parameters_local(15, j_scan) = q.B20_residual;
    parameters_local(16, j_scan) = q.standard_deviation_of_R;
    parameters_local(17, j_scan) = q.standard_deviation_of_Z;
    parameters_local(18, j_scan) = q.beta;

    int_parameters_local[0 + n_int_parameters * j_scan] = q.helicity;
    
    for (j = 0; j < axis_nmax_plus_1; j++) {
      fourier_parameters_local(j + 0 * axis_nmax_plus_1, j_scan) = q.R0c[j];
      fourier_parameters_local(j + 1 * axis_nmax_plus_1, j_scan) = q.R0s[j];
      fourier_parameters_local(j + 2 * axis_nmax_plus_1, j_scan) = q.Z0c[j];
      fourier_parameters_local(j + 3 * axis_nmax_plus_1, j_scan) = q.Z0s[j];
    }
    
    j_scan++;

    if (j_scan >= max_keep_per_proc) keep_going = false;
  }

  end_time = std::chrono::steady_clock::now();
  elapsed = end_time - start_time;
  std::cout << "Proc " << mpi_rank << " finished after " << elapsed.count() << " seconds" << std::endl;

  collect_results(n_parameters, parameters_local, fourier_parameters_local,
		  n_int_parameters, int_parameters_local, j_scan);

}
