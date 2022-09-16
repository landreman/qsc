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
	q.resize_axis_arrays(nf, 0.0);
	
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
	  q.init();
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

///////////////////////////////////////////////////
///////////////////////////////////////////////////

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

  scan1.min_DMerc_times_r2_to_keep = 0;
  scan2.min_DMerc_times_r2_to_keep = scan1.min_DMerc_times_r2_to_keep;
  
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
	for (j = 0; j < qsc::N_FILTERS; j++) {
	  CAPTURE(j);
	  CHECK(scan1.filters[j] == scan2.filters[j]);
	}
	
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
///////////////////////////////////////////////////
///////////////////////////////////////////////////

TEST_CASE("Verify scan results with and without filters are related as expected. [mpi]") {
  int j, k, n;
  qsc::qscfloat amplitude;
  int mpi_rank, n_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  bool proc0 = (mpi_rank == 0);

  // Set up two scans. scan1 will have minimal filters, scan 2 will have an extra filter.
  qsc::Scan scan1;
  qsc::Scan scan2;
  scan1.keep_all = true;
  scan2.keep_all = false;

  scan1.deterministic = true;
  scan2.deterministic = scan1.deterministic;
  
  scan2.max_attempts_per_proc = 60 / n_procs; // Note integer division
  scan1.max_attempts_per_proc = scan2.max_attempts_per_proc;
  
  scan1.max_keep_per_proc = 1000;
  scan2.max_keep_per_proc = scan1.max_keep_per_proc;
  
  scan1.max_seconds = 30;
  scan2.max_seconds = scan1.max_seconds;
  
  scan1.max_d2_volume_d_psi2_to_keep = 0;
  scan2.max_d2_volume_d_psi2_to_keep = scan1.max_d2_volume_d_psi2_to_keep;

  scan1.min_DMerc_times_r2_to_keep = 0;
  scan2.min_DMerc_times_r2_to_keep = scan1.min_DMerc_times_r2_to_keep;
  
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
            
  // Try both O(r^1) and O(r^2):
  for (int order = 1; order < 3; order++) {
    CAPTURE(order);
    if (order == 1) {
      scan1.q.order_r_option = "r1";
    } else if (order == 2) {
      scan1.q.order_r_option = "r2";
    } else {
      throw std::runtime_error("Should not get here");
    }
    scan2.q.order_r_option = scan1.q.order_r_option;
    
    if (order == 2) {
      // Stay close to axisymmetry so O(r^2) diagnostics are not crazy
      scan1.R0c_min[0] = 0.8;
      scan1.R0c_max[0] = 1.2;
      amplitude = 0.02;
    } else {
      // O(r^1). Try some crazier cases
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

    // Run the scan with no filters:
    if (proc0) std::cout << "#### Running scan with keep_all" << std::endl;
    scan1.random();

    int n_filters = 9;
    for (int j_filter = 0; j_filter <= n_filters; j_filter++) {
      CAPTURE(j_filter);

      if (j_filter == 1) {
	scan2.min_R0_to_keep = 0.1;
      } else {
	scan2.min_R0_to_keep = -1.0e+30;
      }

      if (j_filter == 2) {
	scan2.min_iota_to_keep = 0.04;
      } else {
	scan2.min_iota_to_keep = -1.0;
      }

      if (j_filter == 3) {
	scan2.max_elongation_to_keep = 50;
      } else {
	scan2.max_elongation_to_keep = 1.0e+30;
      }

      if (j_filter == 4) {
	scan2.min_L_grad_B_to_keep = 0.25;
      } else {
	scan2.min_L_grad_B_to_keep = -1;
      }

      if (j_filter == 5) {
	scan2.min_L_grad_grad_B_to_keep = 0.5;
      } else {
	scan2.min_L_grad_grad_B_to_keep = -1;
      }

      if (j_filter == 6) {
	scan2.max_B20_variation_to_keep = 0.7;
      } else {
	scan2.max_B20_variation_to_keep = 1.0e+30;
      }

      if (j_filter == 7) {
	scan2.min_r_singularity_to_keep = 0.05;
      } else {
	scan2.min_r_singularity_to_keep = -1;
      }

      if (j_filter == 8) {
	scan2.max_d2_volume_d_psi2_to_keep = 0;
      } else {
	scan2.max_d2_volume_d_psi2_to_keep = 1.0e+30;
      }

      if (j_filter == 9) {
	scan2.min_DMerc_times_r2_to_keep = 0;
      } else { 
	scan2.min_DMerc_times_r2_to_keep = -1.0e+30;
      }
      	
      // Run the scan
      if (proc0) std::cout << "#### Running scan with order=" << order
			   << ", j_filter=" << j_filter << std::endl;
      scan2.random();
      
      // Restore printing format:
      std::cout << std::setprecision(15);
      
      if (proc0) {
	n = -1;
	for (j = 0; j < scan1.n_scan; j++) {
	  CAPTURE(j);
	  CAPTURE(n);
	  // If O(r^1) filters are satisfied:
	  if (scan1.scan_min_R0[j] >= scan2.min_R0_to_keep &&
	      std::abs(scan1.scan_iota[j]) >= scan2.min_iota_to_keep &&
	      scan1.scan_max_elongation[j] <= scan2.max_elongation_to_keep &&
	      scan1.scan_min_L_grad_B[j] >= scan2.min_L_grad_B_to_keep &&
	      (order == 1 || // If O(r^2) filters are satisfied:
	       (scan1.scan_min_L_grad_grad_B[j] >= scan2.min_L_grad_grad_B_to_keep &&
		scan1.scan_B20_variation[j] <= scan2.max_B20_variation_to_keep &&
		scan1.scan_r_singularity[j] >= scan2.min_r_singularity_to_keep &&
		scan1.scan_d2_volume_d_psi2[j] <= scan2.max_d2_volume_d_psi2_to_keep &&
		scan1.scan_DMerc_times_r2[j] >= scan2.min_DMerc_times_r2_to_keep))) {
	    
	    n++;
	    CHECK(Approx(scan1.scan_eta_bar[j]) == scan2.scan_eta_bar[n]);
	    CHECK(Approx(scan1.scan_sigma0[j]) == scan2.scan_sigma0[n]);
	    CHECK(Approx(scan1.scan_B2s[j]) == scan2.scan_B2s[n]);
	    CHECK(Approx(scan1.scan_B2c[j]) == scan2.scan_B2c[n]);
	    CHECK(Approx(scan1.scan_min_R0[j]) == scan2.scan_min_R0[n]);
	    CHECK(Approx(scan1.scan_max_curvature[j]) == scan2.scan_max_curvature[n]);
	    CHECK(Approx(scan1.scan_iota[j]) == scan2.scan_iota[n]);
	    CHECK(Approx(scan1.scan_max_elongation[j]) == scan2.scan_max_elongation[n]);
	    CHECK(Approx(scan1.scan_min_L_grad_B[j]) == scan2.scan_min_L_grad_B[n]);
	    CHECK(Approx(scan1.scan_min_L_grad_grad_B[j]) == scan2.scan_min_L_grad_grad_B[n]);
	    CHECK(Approx(scan1.scan_r_singularity[j]) == scan2.scan_r_singularity[n]);
	    CHECK(Approx(scan1.scan_B20_variation[j]) == scan2.scan_B20_variation[n]);
	    CHECK(Approx(scan1.scan_d2_volume_d_psi2[j]) == scan2.scan_d2_volume_d_psi2[n]);
	    CHECK(Approx(scan1.scan_DMerc_times_r2[j]) == scan2.scan_DMerc_times_r2[n]);
	  
	    CHECK(scan1.scan_helicity[j] == scan2.scan_helicity[n]);
	  
	    for (k = 0; k < nf; k++) {
	      CAPTURE(k);
	      CHECK(Approx(scan1.scan_R0c(k, j)) == scan2.scan_R0c(k, n));
	      CHECK(Approx(scan1.scan_R0s(k, j)) == scan2.scan_R0s(k, n));
	      CHECK(Approx(scan1.scan_Z0c(k, j)) == scan2.scan_Z0c(k, n));
	      CHECK(Approx(scan1.scan_Z0s(k, j)) == scan2.scan_Z0s(k, n));
	    }
	  }
	}
      }
    }
  }
}

///////////////////////////////////////////////////
///////////////////////////////////////////////////

TEST_CASE("Verify the number of configurations attempted or kept in a scan matches the request. [mpi]") {
  int j, k;
  qsc::qscfloat amplitude;
  int mpi_rank, n_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  bool proc0 = (mpi_rank == 0);

  qsc::Scan scan;

  scan.deterministic = true;
    
  scan.max_seconds = 30;
  
  scan.q.nfp = 3;  
  scan.q.nphi = 31;
  scan.q.verbose = 0;  
  scan.q.p2 = -1.0e+4; // Include nonzero pressure so DMerc is nonzero.
  
  int nf = 2;
  scan.R0c_min.resize(nf, 0.0);
  scan.R0c_max.resize(nf, 0.0);
  scan.R0s_min.resize(nf, 0.0);
  scan.R0s_max.resize(nf, 0.0);
  scan.Z0c_min.resize(nf, 0.0);
  scan.Z0c_max.resize(nf, 0.0);
  scan.Z0s_min.resize(nf, 0.0);
  scan.Z0s_max.resize(nf, 0.0);
  
  // Try both O(r^1) and O(r^2):
  for (int order = 1; order < 3; order++) {
    CAPTURE(order);
    // Try both keep_all = true and false:
    for (int j_keep_all = 0; j_keep_all < 2; j_keep_all++) {
      CAPTURE(j_keep_all);
      if (order == 1) {
	scan.q.order_r_option = "r1";
      } else if (order == 2) {
	scan.q.order_r_option = "r2";
      } else {
	throw std::runtime_error("Should not get here");
      }
      scan.q.order_r_option = scan.q.order_r_option;
      
      amplitude = 0.02;
      
      scan.R0c_min[0] = 0.8;
      scan.R0c_max[0] = 1.2;
      
      scan.R0c_min[1] = -amplitude;
      scan.R0c_max[1] =  amplitude;
      
      scan.R0s_min[1] = -amplitude;
      scan.R0s_max[1] =  amplitude;
      
      scan.Z0c_min[1] = -amplitude;
      scan.Z0c_max[1] =  amplitude;
      
      scan.Z0s_min[1] = -amplitude;
      scan.Z0s_max[1] =  amplitude;
      
      scan.eta_bar_min = 0.7;
      scan.eta_bar_max = 1.4;
      
      scan.sigma0_min = -0.3;
      scan.sigma0_max = 0.6;
      
      scan.B2c_min = -0.3;
      scan.B2c_max = 0.3;
      
      scan.B2s_min = -0.3;
      scan.B2s_max = 0.3;
            
      scan.keep_all = (bool) j_keep_all;

      for (int j_constraint = 0; j_constraint < 2; j_constraint++) {
	CAPTURE(j_constraint);
	if (j_constraint == 0) {
	  // Constraints are so lax that they should not be active
	  scan.max_elongation_to_keep = 1.0e+30;
	  scan.min_R0_to_keep = 0.01;
	  scan.max_d2_volume_d_psi2_to_keep = 1.0e+30;
	  scan.min_DMerc_times_r2_to_keep = -1.0e+30;
	  scan.min_L_grad_grad_B_to_keep = -1.0;
	} else {
	  // Constraints should be active
	  scan.max_elongation_to_keep = 2.8;
	  scan.min_R0_to_keep = 0.81;
	  scan.max_d2_volume_d_psi2_to_keep = 0;
	  scan.min_DMerc_times_r2_to_keep = 0;  
	  scan.min_L_grad_grad_B_to_keep = 0.03;
	}
	for (int j_limit = 0; j_limit < 2; j_limit++) {
	  CAPTURE(j_limit);
	  if (j_limit == 0) {
	    // Limited by max_attempts_per_proc
	    scan.max_attempts_per_proc = 30;  
	    scan.max_keep_per_proc = 1000;
	  } else {
	    // Limited by max_keep_per_proc
	    scan.max_attempts_per_proc = 10000;
	    scan.max_keep_per_proc = 4;
	  }

	  for (int j_deterministic = 0; j_deterministic < 2; j_deterministic++) {
	    CAPTURE(j_deterministic);
	    scan.deterministic = (bool) j_deterministic;
	    
	    // Run the scans (without reading an input file):
	    scan.random();
	  
	    if (proc0) {
	      if (j_limit == 0) {
		// Limited by max_attempts_per_proc
		CHECK(scan.filters[qsc::ATTEMPTS] == scan.max_attempts_per_proc * n_procs);
		if (scan.keep_all) CHECK(scan.n_scan == scan.max_attempts_per_proc * n_procs);
	      } else {
		// Limited by max_keep_per_proc
		CHECK(scan.n_scan == scan.max_keep_per_proc * n_procs);
	      }
	    }
	  }
	}
      }
    }
  }
}
