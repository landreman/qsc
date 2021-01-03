#include <chrono>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "qsc.hpp"
#include "scan.hpp"
#include "random.hpp"

using namespace qsc;

/** Given results on all the MPI processes, send them to proc 0.
 */
void Scan::collect_results(int n_parameters,
			   Matrix& parameters_local,
			   Matrix& fourier_parameters_local,
			   int n_int_parameters,
			   std::valarray<int>& int_parameters_local,
			   big j_scan) {
  const int axis_nmax_plus_1 = R0c_max.size();
  const int n_fourier_parameters = axis_nmax_plus_1 * 4;
  int j, k;
  int mpi_rank, n_procs;
  MPI_Status mpi_status;
  
  // Get MPI info
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Comm_size(mpi_comm, &n_procs);
  bool proc0 = (mpi_rank == 0);

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
    std::cout << "Transferring results to proc 0" << std::endl;

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
    std::cout << "Total # of attempts: " << filters[ATTEMPTS] << std::endl;
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

    big total_rejected = 0;
    for (j = KEPT + 1; j < N_FILTERS; j++) total_rejected += filters[j];

    std::cout << "Summary of scan results:" << std::endl;
    std::cout << "  Configurations attempted:          " << filters[ATTEMPTS] << std::endl;
    std::cout << "  # solves of the sigma equation:    " << filters[N_SIGMA_EQ_SOLVES] << std::endl;
    std::cout << "  # solves of the O(r^2) equations:  " << filters[N_R2_SOLVES] << std::endl;
    std::cout << "  Rejected due to crude R0 check:    " << filters[REJECTED_DUE_TO_R0_CRUDE] << std::endl;
    std::cout << "  Rejected due to min_R0:            " << filters[REJECTED_DUE_TO_R0] << std::endl;
    std::cout << "  Rejected due to max curvature:     " << filters[REJECTED_DUE_TO_CURVATURE] << std::endl;
    std::cout << "  Rejected due to min iota:          " << filters[REJECTED_DUE_TO_IOTA] << std::endl;
    std::cout << "  Rejected due to max elongation:    " << filters[REJECTED_DUE_TO_ELONGATION] << std::endl;
    std::cout << "  Rejected due to min L_grad_B:      " << filters[REJECTED_DUE_TO_L_GRAD_B] << std::endl;
    std::cout << "  Rejected due to B20 variation:     " << filters[REJECTED_DUE_TO_B20_VARIATION] << std::endl;
    std::cout << "  Rejected due to min L_grad_grad_B: " << filters[REJECTED_DUE_TO_L_GRAD_GRAD_B] << std::endl;
    std::cout << "  Rejected due to d2_volume_d_psi2:  " << filters[REJECTED_DUE_TO_D2_VOLUME_D_PSI2] << std::endl;
    std::cout << "  Rejected due to DMerc:             " << filters[REJECTED_DUE_TO_DMERC] << std::endl;
    std::cout << "  Rejected due to r_singularity:     " << filters[REJECTED_DUE_TO_R_SINGULARITY] << std::endl;
    std::cout << "  Total rejected:                    " << total_rejected << std::endl;
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
}
