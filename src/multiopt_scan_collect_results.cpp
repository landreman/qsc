#include <chrono>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "multiopt_scan.hpp"

using namespace qsc;

/** Given results on all the MPI processes, send them to proc 0.
 */
void MultiOptScan::collect_results(int n_parameters,
			   Matrix& parameters_local,
			   Matrix& fourier_parameters_local,
			   int n_int_parameters,
			   std::valarray<int>& int_parameters_local,
			   big j_scan) {
  const int n_fourier_parameters = axis_nmax_plus_1 * 4;
  int j, k;
  MPI_Status mpi_status;
  
  auto start = std::chrono::steady_clock::now();

  filters_local[KEPT] = j_scan;

  MPI_Barrier(mpi_comm);

  if (!proc0) {
    // Procs other than 0: Send all results to proc 0

    // Note mpi_rank is used as the tag
    MPI_Send(           &filters_local[0],                     N_FILTERS, MPI_UNSIGNED_LONG_LONG, 0, mpi_rank, mpi_comm);
    MPI_Send(        &parameters_local[0],         n_parameters * j_scan,           MPI_QSCFLOAT, 0, mpi_rank, mpi_comm);
    MPI_Send(&fourier_parameters_local[0], n_fourier_parameters * j_scan,           MPI_QSCFLOAT, 0, mpi_rank, mpi_comm);
    MPI_Send(    &int_parameters_local[0],     n_int_parameters * j_scan,                MPI_INT, 0, mpi_rank, mpi_comm);
    
  } else {
    // Proc 0 does this following block.
    
    std::cout << "##################################################################################" << std::endl;
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Transferring results to proc 0 at time " << elapsed.count() << std::endl;

    std::valarray<big> filters_combined(N_FILTERS * n_procs);
    std::valarray<big> n_solves_kept(n_procs);
    std::valarray<big> attempts_per_proc(n_procs);
    for (j = 0; j < N_FILTERS; j++) filters_combined[j] = filters_local[j];
    // Receive results from all procs about how many runs were kept, etc
    for (j = 1; j < n_procs; j++) {
      // Use mpi_rank as the tag
      MPI_Recv(&filters_combined[N_FILTERS * j], N_FILTERS,
	       MPI_UNSIGNED_LONG_LONG, j, j, mpi_comm, &mpi_status);
    }
    // Add up numbers from all procs:
    for (j = 0; j < N_FILTERS; j++) filters[j] = 0;
    for (k = 0; k < n_procs; k++) {
      for (j = 0; j < N_FILTERS; j++) {
	filters[j] += filters_combined[j + k * N_FILTERS];
      }
      attempts_per_proc[k] = filters_combined[ATTEMPTS + k * N_FILTERS];
      n_solves_kept[k]     = filters_combined[KEPT + k * N_FILTERS];
    }

    n_scan = filters[KEPT];

    std::cout << "Attempts on each proc:";
    for (j = 0; j < n_procs; j++) std::cout << " " << attempts_per_proc[j];
    std::cout << std::endl;
    std::cout << "# solves kept on each proc:";
    for (j = 0; j < n_procs; j++) std::cout << " " << n_solves_kept[j];
    std::cout << std::endl;

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
	       MPI_QSCFLOAT, j, j, mpi_comm, &mpi_status);
      MPI_Recv(&fourier_parameters(0, offset), n_fourier_parameters * n_solves_kept[j],
	       MPI_QSCFLOAT, j, j, mpi_comm, &mpi_status);
      MPI_Recv(&int_parameters[offset * n_int_parameters], n_int_parameters * n_solves_kept[j],
	       MPI_INT, j, j, mpi_comm, &mpi_status);
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
    scan_B20_residual.resize(n_scan, 0.0);
    scan_standard_deviation_of_R.resize(n_scan, 0.0);
    scan_standard_deviation_of_Z.resize(n_scan, 0.0);

    // Not sure why "* n_int_parameters" is here?
    scan_helicity.resize(n_scan * n_int_parameters, 0);

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
    scan_weight_XY2.resize(n_scan, 0.0);
    scan_weight_XY2Prime.resize(n_scan, 0.0);
    scan_weight_Z2.resize(n_scan, 0.0);
    scan_weight_Z2Prime.resize(n_scan, 0.0);
    scan_weight_XY3.resize(n_scan, 0.0);
    scan_weight_XY3Prime.resize(n_scan, 0.0);
    scan_weight_grad_B.resize(n_scan, 0.0);
    scan_weight_grad_grad_B.resize(n_scan, 0.0);
    scan_weight_r_singularity.resize(n_scan, 0.0);
    scan_weight_axis_length.resize(n_scan, 0.0);
    scan_weight_standard_deviation_of_R.resize(n_scan, 0.0);
    
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
      scan_B20_residual[j]      = parameters(14, j);
      scan_standard_deviation_of_R[j] = parameters(15, j);
      scan_standard_deviation_of_Z[j] = parameters(16, j);

      scan_weight_B20[j]                     = parameters(17, j);
      scan_weight_iota[j]                    = parameters(18, j);
      scan_target_iota[j]                    = parameters(19, j);
      scan_weight_elongation[j]              = parameters(20, j);
      scan_weight_curvature[j]               = parameters(21, j);
      scan_weight_R0[j]                      = parameters(22, j);
      scan_target_min_R0[j]                  = parameters(23, j);
      scan_weight_d2_volume_d_psi2[j]        = parameters(24, j);
      scan_max_d2_volume_d_psi2[j]           = parameters(25, j);
      scan_weight_XY2[j]                     = parameters(26, j);
      scan_weight_XY2Prime[j]                = parameters(27, j);
      scan_weight_Z2[j]                      = parameters(28, j);
      scan_weight_Z2Prime[j]                 = parameters(29, j);
      scan_weight_XY3[j]                     = parameters(30, j);
      scan_weight_XY3Prime[j]                = parameters(31, j);
      scan_weight_grad_B[j]                  = parameters(32, j);
      scan_weight_grad_grad_B[j]             = parameters(33, j);
      scan_weight_r_singularity[j]           = parameters(34, j);
      scan_weight_axis_length[j]             = parameters(35, j);
      scan_weight_standard_deviation_of_R[j] = parameters(36, j);
      
      scan_helicity[j] = int_parameters[0 + j * n_int_parameters];
      
      for (k = 0; k < axis_nmax_plus_1; k++) {
 	scan_R0c(k, j) = fourier_parameters(k + 0 * axis_nmax_plus_1, j);
	scan_R0s(k, j) = fourier_parameters(k + 1 * axis_nmax_plus_1, j);
	scan_Z0c(k, j) = fourier_parameters(k + 2 * axis_nmax_plus_1, j);
	scan_Z0s(k, j) = fourier_parameters(k + 3 * axis_nmax_plus_1, j);
      }
    }

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
      std::cout << std::endl;
      std::cout << "min_L_grad_grad_B: " << scan_min_L_grad_grad_B << std::endl;
      std::cout << std::endl;
      std::cout << "d2_volume_d_psi2: " << scan_d2_volume_d_psi2 << std::endl;
      std::cout << std::endl;
      std::cout << "r_singularity: " << scan_r_singularity << std::endl;
      std::cout << std::endl;
      std::cout << "B20_variation: " << scan_B20_variation << std::endl;
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
    
  } // if proc0

  MPI_Barrier(mpi_comm);
  if (proc0 && verbose) {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for collect_results: "
	      << elapsed.count() << " seconds" << std::endl;
  }
}
