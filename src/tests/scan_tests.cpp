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
