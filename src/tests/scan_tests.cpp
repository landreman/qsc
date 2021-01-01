#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include "doctest.h"
#include "scan.hpp"

using doctest::Approx;

TEST_CASE("Each scan result should match a standalone Qsc. [mpi]") {
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  bool proc0 = (mpi_rank == 0);
  
  // Set up a reasonable scan
  qsc::Scan scan;

  // Try both O(r^1) and O(r^2):
  for (int order = 1; order < 3; order++) {
    CAPTURE(order);
    // Try both keep_all = true and false:
    for (int j_keep_all = 0; j_keep_all < 2; j_keep_all++) {
      CAPTURE(j_keep_all);
      scan.q.nfp = 3;
      scan.q.nphi = 31;
      scan.q.verbose = 0;
      scan.q.p2 = -1.0e+4; // Include nonzero pressure so DMerc is nonzero.
      scan.deterministic = true;
      if (order == 1) {
	scan.q.order_r_option = "r1";
      } else if (order == 2) {
	scan.q.order_r_option = "r2";
      } else {
	throw std::runtime_error("Should not get here");
      }
      
      int nf = 2;
      scan.R0c_min.resize(nf, 0.0);
      scan.R0c_max.resize(nf, 0.0);
      scan.R0s_min.resize(nf, 0.0);
      scan.R0s_max.resize(nf, 0.0);
      scan.Z0c_min.resize(nf, 0.0);
      scan.Z0c_max.resize(nf, 0.0);
      scan.Z0s_min.resize(nf, 0.0);
      scan.Z0s_max.resize(nf, 0.0);
      
      scan.R0c_min[0] = 0.8;
      scan.R0c_max[0] = 1.2;
      
      scan.R0c_min[1] = -0.1;
      scan.R0c_max[1] =  0.1;
      
      scan.R0s_min[1] = -0.1;
      scan.R0s_max[1] =  0.1;
      
      scan.Z0c_min[1] = -0.1;
      scan.Z0c_max[1] =  0.1;
      
      scan.Z0s_min[1] = -0.1;
      scan.Z0s_max[1] =  0.1;
      
      scan.eta_bar_min = 0.7;
      scan.eta_bar_max = 1.4;
      
      scan.sigma0_min = -0.3;
      scan.sigma0_max = 0.6;
      
      scan.B2c_min = -1.0;
      scan.B2c_max = 1.0;
      
      scan.B2s_min = -1.0;
      scan.B2s_max = 1.0;
      
      scan.max_attempts_per_proc = 5;
      scan.max_keep_per_proc = 20;
      scan.max_seconds = 30;
      scan.max_elongation_to_keep = 15;
      scan.keep_all = (bool) j_keep_all;
      
      // Run the scan (without reading an input file):
      scan.random();
      
      // Restore printing format:
      std::cout << std::setprecision(15);
      
      if (proc0) {
	
	// Now set up a standalone Qsc object:
	qsc::Qsc q;
	q.verbose = scan.q.verbose;
	q.nfp = scan.q.nfp;
	q.nphi = scan.q.nphi;
	q.p2 = scan.q.p2;
	q.order_r_option = scan.q.order_r_option;
	q.R0c.resize(nf, 0.0);
	q.R0s.resize(nf, 0.0);
	q.Z0c.resize(nf, 0.0);
	q.Z0s.resize(nf, 0.0);
	
	int j, k;
	for (j = 0; j < scan.n_scan; j++) {
	  CAPTURE(j);
	  // Copy inputs from one of the saved calculations in the scan:
	  q.eta_bar = scan.scan_eta_bar[j];
	  q.sigma0 = scan.scan_sigma0[j];
	  q.B2s = scan.scan_B2s[j];
	  q.B2c = scan.scan_B2c[j];
	  for (k = 0; k < nf; k++) {
	    q.R0c[k] = scan.scan_R0c(k, j);
	    q.R0s[k] = scan.scan_R0s(k, j);
	    q.Z0c[k] = scan.scan_Z0c(k, j);
	    q.Z0s[k] = scan.scan_Z0s(k, j);
	  }
	  
	  // Run the standalone calculation
	  q.calculate();
	  
	  CHECK(Approx(q.grid_min_R0) == scan.scan_min_R0[j]);
	  CHECK(Approx(q.grid_max_curvature) == scan.scan_max_curvature[j]);
	  CHECK(Approx(q.iota) == scan.scan_iota[j]);
	  CHECK(Approx(q.grid_max_elongation) == scan.scan_max_elongation[j]);
	  CHECK(Approx(q.grid_min_L_grad_B) == scan.scan_min_L_grad_B[j]);
	  CHECK(Approx(q.grid_min_L_grad_grad_B) == scan.scan_min_L_grad_grad_B[j]);
	  CHECK(Approx(q.r_singularity_robust) == scan.scan_r_singularity[j]);
	  CHECK(Approx(q.d2_volume_d_psi2) == scan.scan_d2_volume_d_psi2[j]);
	  CHECK(Approx(q.DMerc_times_r2) == scan.scan_DMerc_times_r2[j]);
	  CHECK(Approx(q.B20_grid_variation) == scan.scan_B20_variation[j]);
	}
      }
    }
  }
}

TEST_CASE("Verify results of a deterministic scan are independent of number of mpi procs. [mpi]") {
  int j, k;
  qsc::qscfloat amplitude;
  int mpi_rank, n_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  bool proc0 = (mpi_rank == 0);

  // Set up two scans. scan1 will be on 1 proc, scan2 will be on all procs.
  qsc::Scan scan1;
  qsc::Scan scan2;
  scan1.mpi_comm = MPI_COMM_SELF;

  scan1.deterministic = true;
  scan2.deterministic = scan1.deterministic;
  
  // Scans must be limited by max_attempts_per_proc rather than by time or max_keep
  scan2.max_attempts_per_proc = 30;
  scan1.max_attempts_per_proc = scan2.max_attempts_per_proc * n_procs;
  
  scan1.max_keep_per_proc = 1000;
  scan2.max_keep_per_proc = scan1.max_keep_per_proc;
  
  scan1.max_seconds = 30;
  scan2.max_seconds = scan1.max_seconds;
  
  scan1.max_elongation_to_keep = 100;
  scan2.max_elongation_to_keep = scan1.max_elongation_to_keep;

  scan1.min_R0_to_keep = 0.1;
  scan2.min_R0_to_keep = scan1.min_R0_to_keep;

  scan1.max_d2_volume_d_psi2_to_keep = 0;
  scan2.max_d2_volume_d_psi2_to_keep = scan1.max_d2_volume_d_psi2_to_keep;

  scan1.min_DMerc_to_keep = 0;
  scan2.min_DMerc_to_keep = scan1.min_DMerc_to_keep;
  
  scan1.min_L_grad_grad_B_to_keep = 0.05;
  scan2.min_L_grad_grad_B_to_keep = scan1.min_L_grad_grad_B_to_keep;
  
  scan1.q.nfp = 3;
  scan2.q.nfp = scan1.q.nfp;
  
  scan1.q.nphi = 31;
  scan2.q.nphi = scan1.q.nphi;
  
  scan1.q.verbose = 0;
  scan2.q.verbose = scan1.q.verbose;
  
  scan1.q.p2 = -1.0e+4; // Include nonzero pressure so DMerc is nonzero.
  scan2.q.p2 = scan1.q.p2;
  
  int nf = 2;
  scan1.R0c_min.resize(nf, 0.0);
  scan1.R0c_max.resize(nf, 0.0);
  scan1.R0s_min.resize(nf, 0.0);
  scan1.R0s_max.resize(nf, 0.0);
  scan1.Z0c_min.resize(nf, 0.0);
  scan1.Z0c_max.resize(nf, 0.0);
  scan1.Z0s_min.resize(nf, 0.0);
  scan1.Z0s_max.resize(nf, 0.0);
  
  scan2.R0c_min.resize(nf, 0.0);
  scan2.R0c_max.resize(nf, 0.0);
  scan2.R0s_min.resize(nf, 0.0);
  scan2.R0s_max.resize(nf, 0.0);
  scan2.Z0c_min.resize(nf, 0.0);
  scan2.Z0c_max.resize(nf, 0.0);
  scan2.Z0s_min.resize(nf, 0.0);
  scan2.Z0s_max.resize(nf, 0.0);
  
  // Try both O(r^1) and O(r^2):
  for (int order = 1; order < 3; order++) {
    CAPTURE(order);
    // Try both keep_all = true and false:
    for (int j_keep_all = 0; j_keep_all < 2; j_keep_all++) {
      CAPTURE(j_keep_all);
      if (order == 1) {
	scan1.q.order_r_option = "r1";
      } else if (order == 2) {
	scan1.q.order_r_option = "r2";
      } else {
	throw std::runtime_error("Should not get here");
      }
      scan2.q.order_r_option = scan1.q.order_r_option;
      
      if ((bool)j_keep_all) {
	scan1.R0c_min[0] = 0.8;
	scan1.R0c_max[0] = 1.2;
	amplitude = 0.02;
      } else {
	// Try some crazier cases
	scan1.R0c_min[0] = 0.6;
	scan1.R0c_max[0] = 1.2;
	amplitude = 0.8;
      }
      
      scan1.R0c_min[1] = -amplitude;
      scan1.R0c_max[1] =  amplitude;
      
      scan1.R0s_min[1] = -amplitude;
      scan1.R0s_max[1] =  amplitude;
      
      scan1.Z0c_min[1] = -amplitude;
      scan1.Z0c_max[1] =  amplitude;
      
      scan1.Z0s_min[1] = -amplitude;
      scan1.Z0s_max[1] =  amplitude;
      
      scan2.R0c_min = scan1.R0c_min;
      scan2.R0c_max = scan1.R0c_max;
      scan2.R0s_min = scan1.R0s_min;
      scan2.R0s_max = scan1.R0s_max;
      scan2.Z0c_min = scan1.Z0c_min;
      scan2.Z0c_max = scan1.Z0c_max;
      scan2.Z0s_min = scan1.Z0s_min;
      scan2.Z0s_max = scan1.Z0s_max;
      
      scan1.eta_bar_min = 0.7;
      scan1.eta_bar_max = 1.4;
      scan2.eta_bar_min = scan1.eta_bar_min;
      scan2.eta_bar_max = scan1.eta_bar_max;
      
      scan1.sigma0_min = -0.3;
      scan1.sigma0_max = 0.6;
      scan2.sigma0_min = scan1.sigma0_min;
      scan2.sigma0_max = scan1.sigma0_max;
      
      scan1.B2c_min = -1.0;
      scan1.B2c_max = 1.0;
      scan2.B2c_min = scan1.B2c_min;
      scan2.B2c_max = scan1.B2c_max;
      
      scan1.B2s_min = -1.0;
      scan1.B2s_max = 1.0;
      scan2.B2s_min = scan1.B2s_min;
      scan2.B2s_max = scan1.B2s_max;
            
      scan1.keep_all = (bool) j_keep_all;
      scan2.keep_all = scan1.keep_all;
      
      // Run the scans (without reading an input file):
      if (proc0) scan1.random();
      
      MPI_Barrier(MPI_COMM_WORLD); // This might help output look nicer?
      scan2.random();
      
      // Restore printing format:
      std::cout << std::setprecision(15);
      
      if (proc0) {
	CHECK(scan1.n_scan == scan2.n_scan);
	CHECK(scan1.attempts == scan2.attempts);
	CHECK(scan1.rejected_due_to_R0_crude == scan2.rejected_due_to_R0_crude);
	CHECK(scan1.rejected_due_to_R0 == scan2.rejected_due_to_R0);
	CHECK(scan1.rejected_due_to_curvature == scan2.rejected_due_to_curvature);
	CHECK(scan1.rejected_due_to_iota == scan2.rejected_due_to_iota);
	CHECK(scan1.rejected_due_to_elongation == scan2.rejected_due_to_elongation);
	CHECK(scan1.rejected_due_to_L_grad_B == scan2.rejected_due_to_L_grad_B);
	CHECK(scan1.rejected_due_to_L_grad_grad_B == scan2.rejected_due_to_L_grad_grad_B);
	CHECK(scan1.rejected_due_to_B20_variation == scan2.rejected_due_to_B20_variation);
	CHECK(scan1.rejected_due_to_r_singularity == scan2.rejected_due_to_r_singularity);
	CHECK(scan1.rejected_due_to_d2_volume_d_psi2 == scan2.rejected_due_to_d2_volume_d_psi2);
	CHECK(scan1.rejected_due_to_DMerc == scan2.rejected_due_to_DMerc);
	
	for (j = 0; j < scan1.n_scan; j++) {
	  CAPTURE(j);
	  CHECK(Approx(scan1.scan_eta_bar[j]) == scan2.scan_eta_bar[j]);
	  CHECK(Approx(scan1.scan_sigma0[j]) == scan2.scan_sigma0[j]);
	  CHECK(Approx(scan1.scan_B2s[j]) == scan2.scan_B2s[j]);
	  CHECK(Approx(scan1.scan_B2c[j]) == scan2.scan_B2c[j]);
	  CHECK(Approx(scan1.scan_min_R0[j]) == scan2.scan_min_R0[j]);
	  CHECK(Approx(scan1.scan_max_curvature[j]) == scan2.scan_max_curvature[j]);
	  CHECK(Approx(scan1.scan_iota[j]) == scan2.scan_iota[j]);
	  CHECK(Approx(scan1.scan_max_elongation[j]) == scan2.scan_max_elongation[j]);
	  CHECK(Approx(scan1.scan_min_L_grad_B[j]) == scan2.scan_min_L_grad_B[j]);
	  CHECK(Approx(scan1.scan_min_L_grad_grad_B[j]) == scan2.scan_min_L_grad_grad_B[j]);
	  CHECK(Approx(scan1.scan_r_singularity[j]) == scan2.scan_r_singularity[j]);
	  CHECK(Approx(scan1.scan_B20_variation[j]) == scan2.scan_B20_variation[j]);
	  CHECK(Approx(scan1.scan_d2_volume_d_psi2[j]) == scan2.scan_d2_volume_d_psi2[j]);
	  CHECK(Approx(scan1.scan_DMerc_times_r2[j]) == scan2.scan_DMerc_times_r2[j]);
	  
	  CHECK(scan1.scan_helicity[j] == scan2.scan_helicity[j]);
	  
	  for (k = 0; k < nf; k++) {
	    CAPTURE(k);
	    CHECK(Approx(scan1.scan_R0c(k, j)) == scan2.scan_R0c(k, j));
	    CHECK(Approx(scan1.scan_R0s(k, j)) == scan2.scan_R0s(k, j));
	    CHECK(Approx(scan1.scan_Z0c(k, j)) == scan2.scan_Z0c(k, j));
	    CHECK(Approx(scan1.scan_Z0s(k, j)) == scan2.scan_Z0s(k, j));
	  }
	}
      }
    }
  }
}
