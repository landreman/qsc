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
  const int n_parameters = 14;
  const int n_int_parameters = 1;
  const int axis_nmax_plus_1 = R0c_max.size();
  const int n_big_parameters = 13;
  const int n_fourier_parameters = axis_nmax_plus_1 * 4;
  Matrix parameters_local(n_parameters, max_keep_per_proc);
  Matrix fourier_parameters_local(n_fourier_parameters, max_keep_per_proc);
  std::valarray<int> int_parameters_local(n_int_parameters * max_keep_per_proc);
  std::valarray<big> big_parameters_local(n_big_parameters), big_parameters;
  big j_scan = 0, attempts_local = 0;
  int j, k;
  qscfloat R0_at_0, R0_at_half_period, val;
  int mpi_rank, n_procs;
  MPI_Comm MPI_COMM_QSC = MPI_COMM_WORLD;
  MPI_Status mpi_status;
  
  // Initialize MPI
  MPI_Comm_rank(MPI_COMM_QSC, &mpi_rank);
  MPI_Comm_size(MPI_COMM_QSC, &n_procs);
  bool proc0 = (mpi_rank == 0);
  
  // Initialize random distributions:
  Random random_eta_bar(deterministic, eta_bar_scan_option, eta_bar_min, eta_bar_max);
  Random random_sigma0(deterministic, sigma0_scan_option, sigma0_min, sigma0_max);
  Random random_B2c(deterministic, B2c_scan_option, B2c_min, B2c_max);
  Random random_B2s(deterministic, B2s_scan_option, B2s_min, B2s_max);
  random_eta_bar.set_to_nth(mpi_rank * max_attempts_per_proc);
  random_sigma0.set_to_nth(mpi_rank * max_attempts_per_proc);
  random_B2c.set_to_nth(mpi_rank * max_attempts_per_proc);
  random_B2s.set_to_nth(mpi_rank * max_attempts_per_proc);
  Random* random_R0c[axis_nmax_plus_1];
  Random* random_R0s[axis_nmax_plus_1];
  Random* random_Z0c[axis_nmax_plus_1];
  Random* random_Z0s[axis_nmax_plus_1];
  for (j = 0; j < axis_nmax_plus_1; j++) {
    random_R0c[j] = new Random(deterministic, fourier_scan_option, R0c_min[j], R0c_max[j]);
    random_R0s[j] = new Random(deterministic, fourier_scan_option, R0s_min[j], R0s_max[j]);
    random_Z0c[j] = new Random(deterministic, fourier_scan_option, Z0c_min[j], Z0c_max[j]);
    random_Z0s[j] = new Random(deterministic, fourier_scan_option, Z0s_min[j], Z0s_max[j]);
    random_R0c[j]->set_to_nth(mpi_rank * max_attempts_per_proc);
    random_R0s[j]->set_to_nth(mpi_rank * max_attempts_per_proc + 1);
    random_Z0c[j]->set_to_nth(mpi_rank * max_attempts_per_proc + 2);
    random_Z0s[j]->set_to_nth(mpi_rank * max_attempts_per_proc + 3);
    // Above, the +1 .. +3 is so R0c/R0s/Z0c/Z0s are a bit less correlated
  }
  
  rejected_due_to_R0_crude = 0;
  rejected_due_to_R0 = 0;
  rejected_due_to_curvature = 0;
  rejected_due_to_iota = 0;
  rejected_due_to_elongation = 0;
  rejected_due_to_L_grad_B = 0;
  rejected_due_to_L_grad_grad_B = 0;
  rejected_due_to_B20_variation = 0;
  rejected_due_to_r_singularity = 0;
  rejected_due_to_d2_volume_d_psi2 = 0;
  rejected_due_to_DMerc = 0;

  q.R0c.resize(axis_nmax_plus_1, 0.0);
  q.R0s.resize(axis_nmax_plus_1, 0.0);
  q.Z0c.resize(axis_nmax_plus_1, 0.0);
  q.Z0s.resize(axis_nmax_plus_1, 0.0);

  std::chrono::time_point<std::chrono::steady_clock> start_time, end_time;
  start_time = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed;

  // Initialize the Qsc object:
  q.validate();
  q.allocate();
  
  bool keep_going = true;
  while (keep_going) {
    if (attempts_local >= max_attempts_per_proc) break;
    
    attempts_local++;

    // Initialize R0c, and do a crude check of whether R0 goes negative:
    R0_at_half_period = 0;
    for (j = 0; j < axis_nmax_plus_1; j++) {
      val = random_R0c[j]->get();
      q.R0c[j] = val;
      if (j % 2 == 0) {
	R0_at_half_period += val;
      } else {
	R0_at_half_period -= val;
      }
    }
    R0_at_0 = q.R0c.sum();
    if (R0_at_0 <= 0 || R0_at_half_period <= 0) {
      rejected_due_to_R0_crude++;
      continue;
    }

    // Initialize the rest of the axis shape:
    for (j = 0; j < axis_nmax_plus_1; j++) {
      q.R0s[j] = random_R0s[j]->get();
      q.Z0s[j] = random_Z0s[j]->get();
      q.Z0c[j] = random_Z0c[j]->get();
    }
    q.init_axis();
    if (q.grid_min_R0 < min_R0_to_keep) {
      rejected_due_to_R0++;
      continue;
    }
    if (!keep_all && 1.0 / q.grid_max_curvature < min_L_grad_B_to_keep) {
      rejected_due_to_curvature++;
      continue;
    }

    // Initialize the rest of the O(r^1) parameters:
    q.eta_bar = random_eta_bar.get();
    q.sigma0 = random_sigma0.get();

    // Here is the main O(r^1) solve:
    q.solve_sigma_equation();
    q.r1_diagnostics();
    if (!keep_all && q.grid_max_elongation > max_elongation_to_keep) {
      rejected_due_to_elongation++;
      continue;
    }
    if (!keep_all && q.grid_min_L_grad_B < min_L_grad_B_to_keep) {
      rejected_due_to_L_grad_B++;
      continue;
    }

    if (q.at_least_order_r2) {
      // Initialize the rest of the O(r^2) parameters:
      q.B2c = random_B2c.get();
      q.B2s = random_B2s.get();

      // Here is the main O(r^2) solve:
      q.calculate_r2();

      // Filter results:
      if (!keep_all && q.B20_grid_variation > max_B20_variation_to_keep) {
	rejected_due_to_B20_variation++;
	continue;
      }

      q.mercier();
      if (!keep_all && q.d2_volume_d_psi2 > max_d2_volume_d_psi2_to_keep) {
	rejected_due_to_d2_volume_d_psi2++;
	continue;
      }
      if (!keep_all && q.DMerc_times_r2 < min_DMerc_to_keep) {
	rejected_due_to_DMerc++;
	continue;
      }

      q.calculate_grad_grad_B_tensor();
      if (!keep_all && q.grid_min_L_grad_grad_B < min_L_grad_grad_B_to_keep) {
	rejected_due_to_L_grad_grad_B++;
	continue;
      }

      q.calculate_r_singularity();
      if (!keep_all && q.r_singularity_robust < min_r_singularity_to_keep) {
	rejected_due_to_r_singularity++;
	continue;
      }
    } // if at_least_order_r2

    // If we made it this far, then we found a keeper.
    parameters_local(0 , j_scan) = q.eta_bar;
    parameters_local(1 , j_scan) = q.sigma0;
    parameters_local(2 , j_scan) = q.B2c;
    parameters_local(3 , j_scan) = q.B2s;
    parameters_local(4 , j_scan) = q.grid_min_R0;
    parameters_local(5 , j_scan) = q.grid_max_curvature;
    parameters_local(6 , j_scan) = q.iota;
    parameters_local(7 , j_scan) = q.grid_max_elongation;
    parameters_local(8 , j_scan) = q.grid_min_L_grad_B;
    parameters_local(9 , j_scan) = q.grid_min_L_grad_grad_B;
    parameters_local(10, j_scan) = q.r_singularity_robust;
    parameters_local(11, j_scan) = q.d2_volume_d_psi2;
    parameters_local(12, j_scan) = q.DMerc_times_r2;
    parameters_local(13, j_scan) = q.B20_grid_variation;

    int_parameters_local[0 + n_int_parameters * j_scan] = q.helicity;
    
    for (j = 0; j < axis_nmax_plus_1; j++) {
      fourier_parameters_local(j + 0 * axis_nmax_plus_1, j_scan) = q.R0c[j];
      fourier_parameters_local(j + 1 * axis_nmax_plus_1, j_scan) = q.R0s[j];
      fourier_parameters_local(j + 2 * axis_nmax_plus_1, j_scan) = q.Z0c[j];
      fourier_parameters_local(j + 3 * axis_nmax_plus_1, j_scan) = q.Z0s[j];
    }
    
    j_scan++;

    if (j_scan >= max_keep_per_proc) keep_going = false;

    end_time = std::chrono::steady_clock::now();
    elapsed = end_time - start_time;
    if (elapsed.count() > max_seconds) keep_going = false;
  }

  big_parameters_local[ 0] = attempts_local;
  big_parameters_local[ 1] = j_scan;
  big_parameters_local[ 2] = rejected_due_to_R0_crude;
  big_parameters_local[ 3] = rejected_due_to_R0;
  big_parameters_local[ 4] = rejected_due_to_curvature;
  big_parameters_local[ 5] = rejected_due_to_iota;
  big_parameters_local[ 6] = rejected_due_to_elongation;
  big_parameters_local[ 7] = rejected_due_to_L_grad_B;
  big_parameters_local[ 8] = rejected_due_to_B20_variation;
  big_parameters_local[ 9] = rejected_due_to_L_grad_grad_B;
  big_parameters_local[10] = rejected_due_to_d2_volume_d_psi2;
  big_parameters_local[11] = rejected_due_to_DMerc;
  big_parameters_local[12] = rejected_due_to_r_singularity;

  end_time = std::chrono::steady_clock::now();
  elapsed = end_time - start_time;
  std::cout << "Proc " << mpi_rank << " finished after " << elapsed.count() << " seconds" << std::endl;

  MPI_Barrier(MPI_COMM_QSC);

  if (!proc0) {
    // Procs other than 0: Send all results to proc 0

    // Note mpi_rank is used as the tag
    MPI_Send(    &big_parameters_local[0],              n_big_parameters, MPI_UNSIGNED_LONG, 0, mpi_rank, MPI_COMM_QSC);
    MPI_Send(        &parameters_local[0],         n_parameters * j_scan,      MPI_QSCFLOAT, 0, mpi_rank, MPI_COMM_QSC);
    MPI_Send(&fourier_parameters_local[0], n_fourier_parameters * j_scan,      MPI_QSCFLOAT, 0, mpi_rank, MPI_COMM_QSC);
    MPI_Send(    &int_parameters_local[0],     n_int_parameters * j_scan,           MPI_INT, 0, mpi_rank, MPI_COMM_QSC);
    
  } else {
    // Proc 0 does this following block.
    
    std::cout << "##################################################################################" << std::endl;
    std::cout << "Transferring results to proc 0" << std::endl;

    big_parameters.resize(n_big_parameters * n_procs, 0);
    std::valarray<big> n_solves_kept(n_procs);
    std::valarray<big> attempts_per_proc(n_procs);
    for (j = 0; j < n_big_parameters; j++) big_parameters[j] = big_parameters_local[j];
    // Receive results from all procs about how many runs were kept, etc
    for (j = 1; j < n_procs; j++) {
      // Use mpi_rank as the tag
      MPI_Recv(&big_parameters[n_big_parameters * j], n_big_parameters,
	       MPI_UNSIGNED_LONG, j, j, MPI_COMM_QSC, &mpi_status);
    }
    // Add up numbers from all procs:
    std::valarray<big> total_big_parameters;
    total_big_parameters.resize(n_big_parameters, 0);
    for (k = 0; k < n_procs; k++) {
      for (j = 0; j < n_big_parameters; j++) {
	total_big_parameters[j] += big_parameters[j + k * n_big_parameters];
      }
      attempts_per_proc[k] = big_parameters[0 + k * n_big_parameters];
      n_solves_kept[k]     = big_parameters[1 + k * n_big_parameters];
    }

    // Unpack:
    attempts                         = total_big_parameters[0];
    n_scan                           = total_big_parameters[1];
    rejected_due_to_R0_crude         = total_big_parameters[2];
    rejected_due_to_R0               = total_big_parameters[3];
    rejected_due_to_curvature        = total_big_parameters[4];
    rejected_due_to_iota             = total_big_parameters[5];
    rejected_due_to_elongation       = total_big_parameters[6];
    rejected_due_to_L_grad_B         = total_big_parameters[7];
    rejected_due_to_B20_variation    = total_big_parameters[8];
    rejected_due_to_L_grad_grad_B    = total_big_parameters[9];
    rejected_due_to_d2_volume_d_psi2 = total_big_parameters[10];
    rejected_due_to_DMerc            = total_big_parameters[11];
    rejected_due_to_r_singularity    = total_big_parameters[12];

    std::cout << "Attempts on each proc:";
    for (j = 0; j < n_procs; j++) std::cout << " " << attempts_per_proc[j];
    std::cout << std::endl;
    std::cout << "# solves kept on each proc:";
    for (j = 0; j < n_procs; j++) std::cout << " " << n_solves_kept[j];
    std::cout << std::endl;
    std::cout << "Total # of attempts: " << attempts << std::endl;
    std::cout << "Total # of solves kept: " << n_scan << std::endl;

    // Now that we know the total # of runs that succeeded, we can
    // allocate big arrays to store the combined results from all
    // procs.
    Matrix parameters(n_parameters, n_scan);
    Matrix fourier_parameters(n_fourier_parameters, n_scan);
    std::valarray<int> int_parameters(n_int_parameters * n_scan);
    // Copy proc0 results to the final arrays:
    for (j = 0; j < j_scan; j++) {
      for (k = 0; k < n_parameters; k++) {
	parameters(k, j) = parameters_local(k, j);
      }
      for (k = 0; k < n_fourier_parameters; k++) {
	fourier_parameters(k, j) = fourier_parameters_local(k, j);
      }
      for (k = 0; k < n_int_parameters; k++) {
	int_parameters[k + n_int_parameters * j] = int_parameters_local[k + n_int_parameters * j];
      }
    }
    // Receive results from other procs:
    big offset = j_scan;
    for (j = 1; j < n_procs; j++) {
      // Use mpi_rank as the tag
      MPI_Recv(&parameters(0, offset), n_parameters * n_solves_kept[j],
	       MPI_QSCFLOAT, j, j, MPI_COMM_QSC, &mpi_status);
      MPI_Recv(&fourier_parameters(0, offset), n_fourier_parameters * n_solves_kept[j],
	       MPI_QSCFLOAT, j, j, MPI_COMM_QSC, &mpi_status);
      MPI_Recv(&int_parameters[offset * n_int_parameters], n_int_parameters * n_solves_kept[j],
	       MPI_INT, j, j, MPI_COMM_QSC, &mpi_status);
      offset += n_solves_kept[j];
    }

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

    scan_helicity.resize(n_scan * n_int_parameters, 0);

    scan_R0c.resize(axis_nmax_plus_1, n_scan, 0.0);
    scan_R0s.resize(axis_nmax_plus_1, n_scan, 0.0);
    scan_Z0c.resize(axis_nmax_plus_1, n_scan, 0.0);
    scan_Z0s.resize(axis_nmax_plus_1, n_scan, 0.0);
    
    // Unpack parameters
    for (j = 0; j < n_scan; j++) {
      scan_eta_bar[j]           = parameters( 0, j);
      scan_sigma0[j]            = parameters( 1, j);
      scan_B2c[j]               = parameters( 2, j);
      scan_B2s[j]               = parameters( 3, j);
      scan_min_R0[j]            = parameters( 4, j);
      scan_max_curvature[j]     = parameters( 5, j);
      scan_iota[j]              = parameters( 6, j);
      scan_max_elongation[j]    = parameters( 7, j);
      scan_min_L_grad_B[j]      = parameters( 8, j);
      scan_min_L_grad_grad_B[j] = parameters( 9, j);
      scan_r_singularity[j]     = parameters(10, j);
      scan_d2_volume_d_psi2[j]  = parameters(11, j);
      scan_DMerc_times_r2[j]    = parameters(12, j);
      scan_B20_variation[j]     = parameters(13, j);
      
      scan_helicity[j] = int_parameters[0 + j * n_int_parameters];
      
      for (k = 0; k < axis_nmax_plus_1; k++) {
	scan_R0c(k, j) = fourier_parameters(k + 0 * axis_nmax_plus_1, j);
	scan_R0s(k, j) = fourier_parameters(k + 1 * axis_nmax_plus_1, j);
	scan_Z0c(k, j) = fourier_parameters(k + 2 * axis_nmax_plus_1, j);
	scan_Z0s(k, j) = fourier_parameters(k + 3 * axis_nmax_plus_1, j);
      }
    }
  } // end of MPI communication of results.
  
  big total_rejected = rejected_due_to_R0_crude
    + rejected_due_to_R0
    + rejected_due_to_curvature
    + rejected_due_to_iota
    + rejected_due_to_elongation
    + rejected_due_to_L_grad_B
    + rejected_due_to_B20_variation
    + rejected_due_to_L_grad_grad_B
    + rejected_due_to_d2_volume_d_psi2
    + rejected_due_to_DMerc
    + rejected_due_to_r_singularity;

  if (proc0) {
    std::cout << "Summary of scan results:" << std::endl;
    std::cout << "  Configurations attempted:          " << attempts << std::endl;
    std::cout << "  Rejected due to crude R0 check:    " << rejected_due_to_R0_crude << std::endl;
    std::cout << "  Rejected due to min_R0:            " << rejected_due_to_R0 << std::endl;
    std::cout << "  Rejected due to max curvature:     " << rejected_due_to_curvature << std::endl;
    std::cout << "  Rejected due to min iota:          " << rejected_due_to_iota << std::endl;
    std::cout << "  Rejected due to max elongation:    " << rejected_due_to_elongation << std::endl;
    std::cout << "  Rejected due to min L_grad_B:      " << rejected_due_to_L_grad_B << std::endl;
    std::cout << "  Rejected due to B20 variation:     " << rejected_due_to_B20_variation << std::endl;
    std::cout << "  Rejected due to min L_grad_grad_B: " << rejected_due_to_L_grad_grad_B << std::endl;
    std::cout << "  Rejected due to d2_volume_d_psi2:  " << rejected_due_to_d2_volume_d_psi2 << std::endl;
    std::cout << "  Rejected due to DMerc:             " << rejected_due_to_DMerc << std::endl;
    std::cout << "  Rejected due to r_singularity:     " << rejected_due_to_r_singularity << std::endl;
    std::cout << "  Kept:                              " << n_scan << std::endl;
    std::cout << "  Kept + rejected:                   " << n_scan + total_rejected << std::endl;
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
      std::cout << "r_singularity: " << scan_r_singularity << std::endl;
    }
  }

  MPI_Barrier(MPI_COMM_QSC);
}
