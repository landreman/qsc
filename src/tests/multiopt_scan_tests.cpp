#include <iostream>
#include <mpi.h>
#include "doctest.h"
#include "multiopt_scan.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("Check linear and logarithmic spacing of scan parameters. [multiopt_scan]") {
  if (single) return;
  
  MultiOptScan mos;
  mos.mo_ref.opts.resize(1); // mos.init() requires that at least 1 opt be present.
  mos.params = {"R0c1", "weight_B20"};
  mos.params_min = {0.1, 0.01};
  mos.params_max = {0.4, 1.0};
  mos.params_n = {4, 3};
  mos.params_log = {false, true};
  mos.params_stage = {-1, -1};
  
  mos.init();
  CHECK(mos.n_scan_all == 12);
  
  CHECK(Approx(mos.params_vals[0][0]) == 0.1);
  CHECK(Approx(mos.params_vals[0][1]) == 0.2);
  CHECK(Approx(mos.params_vals[0][2]) == 0.3);
  CHECK(Approx(mos.params_vals[0][3]) == 0.4);
  
  CHECK(Approx(mos.params_vals[1][0]) == 0.01);
  CHECK(Approx(mos.params_vals[1][1]) == 0.1);
  CHECK(Approx(mos.params_vals[1][2]) == 1.0);
}

TEST_CASE("Run a small MultiOptScan with keep_all true. [mpi] [multiopt_scan]") {
  // This example runs an optimization for QH with nfp = 4.
  
  if (single) return;

  // This example only works for 2-4 procs.
  int n_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  if (n_procs < 2 || n_procs > 4) return;
      
  MultiOptScan mos;

  // Set scan parameters:
  mos.verbose = 2;
  mos.keep_all = true;
  mos.max_elongation_to_keep = 3.0; // This filter should not be active though since keep_all = true.

  mos.params       = {"R0c1", "weight_B20"};
  mos.params_min   = {  0.15,          0.1};
  mos.params_max   = {   0.2,          1.0};
  mos.params_n     = {     2,            2};
  mos.params_log   = { false,         true};
  mos.params_stage = {    -1,           -1};

  // Set 2 opt stages:
  mos.mo_ref.opts.resize(2);
  mos.mo_ref.verbose = 2;

  // Set the initial QSC configuration:
  mos.mo_ref.opts[0].q.nfp = 4;
  mos.mo_ref.opts[0].q.nphi = 61; // This value should be over-ridden by mo_ref.opts[0].nphi
  mos.mo_ref.opts[0].q.verbose = 0;
  mos.mo_ref.opts[0].q.order_r_option = "r2.1";
  mos.mo_ref.opts[0].q.eta_bar = 1.0;
  mos.mo_ref.opts[0].q.R0c = {1.0, 0.17};
  mos.mo_ref.opts[0].q.Z0s = {0.0, 0.17};
  mos.mo_ref.opts[0].q.R0s = {0.0, 0.0};
  mos.mo_ref.opts[0].q.Z0c = {0.0, 0.0};

  // Set parameters for opt stage 0:
  mos.mo_ref.opts[0].verbose = 0;
  mos.mo_ref.opts[0].fourier_refine = 2;
  mos.mo_ref.opts[0].nphi = {17, 19, 21};
  mos.mo_ref.opts[0].vary_eta_bar = true;
  mos.mo_ref.opts[0].vary_B2c = true;
  mos.mo_ref.opts[0].vary_R0c = {false, false};
  mos.mo_ref.opts[0].vary_Z0s = {false, true};
  mos.mo_ref.opts[0].vary_R0s = {false, false};
  mos.mo_ref.opts[0].vary_Z0c = {false, false};
  mos.mo_ref.opts[0].weight_grad_B = 1.0;
  mos.mo_ref.opts[0].weight_B20 = 1.0;
  
  // Set parameters for opt stage 1:
  mos.mo_ref.opts[1].verbose = 0;
  mos.mo_ref.opts[1].fourier_refine = 0;
  mos.mo_ref.opts[1].nphi = {25};
  mos.mo_ref.opts[1].vary_eta_bar = true;
  mos.mo_ref.opts[1].vary_B2c = true;
  mos.mo_ref.opts[1].vary_R0c = {false, false, true, true};
  mos.mo_ref.opts[1].vary_Z0s = {false, true, true, true};
  mos.mo_ref.opts[1].vary_R0s = {false, false, false, false};
  mos.mo_ref.opts[1].vary_Z0c = {false, false, false, false};
  mos.mo_ref.opts[1].weight_grad_B = 1.0;
  mos.mo_ref.opts[1].weight_B20 = 2.0;
  mos.mo_ref.opts[1].weight_grad_grad_B = 0.01;

  // Run the scan:
  mos.init();
  mos.scan();

  // Check results. Only MPI proc 0 has the final data.
  if (mos.proc0) {
    CHECK(mos.n_scan == 4);

    CHECK(Approx(mos.scan_weight_B20[0]) == 0.1);
    CHECK(Approx(mos.scan_weight_B20[1]) == 1.0);
    CHECK(Approx(mos.scan_weight_B20[2]) == 0.1);
    CHECK(Approx(mos.scan_weight_B20[3]) == 1.0);
    CHECK(Approx(mos.scan_R0c(1, 0)) == 0.15);
    CHECK(Approx(mos.scan_R0c(1, 1)) == 0.15);
    CHECK(Approx(mos.scan_R0c(1, 2)) == 0.2);
    CHECK(Approx(mos.scan_R0c(1, 3)) == 0.2);

    CHECK(Approx(mos.scan_B20_variation[0]) == 2.26940962262729);
    CHECK(Approx(mos.scan_B20_variation[1]) == 0.173028892334032);
    CHECK(Approx(mos.scan_B20_variation[2]) == 0.630219083075497);
    CHECK(Approx(mos.scan_B20_variation[3]) == 0.291004337639454);

  } else {
    // Only procs other than 0 have actually done any solves.
    CHECK(mos.mo.opts[0].q.nphi == 21);
    CHECK(mos.mo.opts[1].q.nphi == 25);
  }

  // The total number of function evals for the scan should match the
  // sum of the number of function evals for each entry in the scan,
  // since all points are kept.
  if (mos.proc0) {
    int total_n_evals = 0;
    for (int j = 0; j < mos.n_scan; j++)
      total_n_evals += mos.scan_n_evals[j];
    CHECK(total_n_evals == mos.n_evals);
  }
  
  // For each entry in the scan, run a standalone single QSC and
  // confirm all the results match.
  if (mos.proc0) {
    for (int j = 0; j < mos.n_scan; j++) {
      CAPTURE(j);
      
      Qsc q;
      q.nfp = mos.mo_ref.opts[0].q.nfp;
      q.nphi = 25; // Must match value above in this test.
      q.verbose = 0;
      q.eta_bar = mos.scan_eta_bar[j];
      q.sigma0 = mos.scan_sigma0[j];
      q.B2c = mos.scan_B2c[j];
      q.B2s = mos.scan_B2s[j];
      q.order_r_option = "r2.1";
      q.resize_axis_arrays(mos.axis_nmax_plus_1, 0.0);
      for (int k = 0; k < mos.axis_nmax_plus_1; k++) {
	q.R0c[k] = mos.scan_R0c(k, j);
	q.R0s[k] = mos.scan_R0s(k, j);
	q.Z0c[k] = mos.scan_Z0c(k, j);
	q.Z0s[k] = mos.scan_Z0s(k, j);
      }

      q.init();
      q.calculate();

      CHECK(Approx(q.grid_min_R0) == mos.scan_min_R0[j]);
      CHECK(Approx(q.grid_max_curvature) == mos.scan_max_curvature[j]);
      CHECK(Approx(q.iota) == mos.scan_iota[j]);
      CHECK(Approx(q.grid_max_elongation) == mos.scan_max_elongation[j]);
      CHECK(Approx(q.grid_min_L_grad_B) == mos.scan_min_L_grad_B[j]);
      CHECK(Approx(q.grid_min_L_grad_grad_B) == mos.scan_min_L_grad_grad_B[j]);
      CHECK(Approx(q.r_singularity_robust) == mos.scan_r_singularity[j]);
      CHECK(Approx(q.B20_grid_variation) == mos.scan_B20_variation[j]);
      CHECK(Approx(q.B20_residual) == mos.scan_B20_residual[j]);
      CHECK(Approx(q.B20_mean) == mos.scan_B20_mean[j]);
      CHECK(Approx(q.d2_volume_d_psi2) == mos.scan_d2_volume_d_psi2[j]);
      CHECK(Approx(q.DMerc_times_r2) == mos.scan_DMerc_times_r2[j]);
      CHECK(Approx(q.standard_deviation_of_R) == mos.scan_standard_deviation_of_R[j]);
      CHECK(Approx(q.standard_deviation_of_Z) == mos.scan_standard_deviation_of_Z[j]);
      CHECK(q.helicity == mos.scan_helicity[j]);
      CHECK(Approx(q.grid_max_XY2) == mos.scan_max_XY2[j]);
      CHECK(Approx(q.grid_max_Z2) == mos.scan_max_Z2[j]);
      CHECK(Approx(q.grid_max_XY3) == mos.scan_max_XY3[j]);
      CHECK(Approx(q.grid_max_d_XY2_d_varphi) == mos.scan_max_d_XY2_d_varphi[j]);
      CHECK(Approx(q.grid_max_d_Z2_d_varphi) == mos.scan_max_d_Z2_d_varphi[j]);
      CHECK(Approx(q.grid_max_d_XY3_d_varphi) == mos.scan_max_d_XY3_d_varphi[j]);
      CHECK(Approx(q.grid_max_d2_XY2_d_varphi2) == mos.scan_max_d2_XY2_d_varphi2[j]);
      CHECK(Approx(q.grid_max_d2_XY3_d_varphi2) == mos.scan_max_d2_XY3_d_varphi2[j]);
    }
  }
}

TEST_CASE("Run a small MultiOptScan with keep_all false. [mpi] [multiopt_scan]") {
  // This example runs an optimization for QA with nfp = 2.
  
  if (single) return;
  
  // This example only works for 2-5 procs.
  int n_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  if (n_procs < 2 || n_procs > 5) return;
  
  MultiOptScan mos;

  // Set scan parameters:
  mos.verbose = 2;
  mos.keep_all = false;
  mos.min_R0_to_keep = 0.5; // This filter should not be active.
  mos.max_elongation_to_keep = 3.0; // This filter should not be active.
  mos.min_r_singularity_to_keep = 0.2;
  mos.max_B20_variation_to_keep = 0.04;

  mos.params       = {"weight_B20", "weight_iota"};
  mos.params_min   = {       0.001,           8.0};
  mos.params_max   = {        12.0,          10.0};
  mos.params_n     = {           2,             2};
  mos.params_log   = {        true,          true};
  mos.params_stage = {           1,            -1};

  // Set 2 opt stages:
  mos.mo_ref.opts.resize(2);
  mos.mo_ref.verbose = 2;

  // Set the initial QSC configuration:
  mos.mo_ref.opts[0].q.nfp = 2;
  mos.mo_ref.opts[0].q.nphi = 25;
  mos.mo_ref.opts[0].q.verbose = 0;
  mos.mo_ref.opts[0].q.order_r_option = "r2.1";
  mos.mo_ref.opts[0].q.eta_bar = 1.0;
  mos.mo_ref.opts[0].q.R0c = {1.0, 0.17};
  mos.mo_ref.opts[0].q.Z0s = {0.0, 0.17};
  mos.mo_ref.opts[0].q.R0s = {0.0, 0.0};
  mos.mo_ref.opts[0].q.Z0c = {0.0, 0.0};

  // Set parameters for opt stage 0:
  mos.mo_ref.opts[0].verbose = 0;
  mos.mo_ref.opts[0].fourier_refine = 2;
  mos.mo_ref.opts[0].vary_eta_bar = true;
  mos.mo_ref.opts[0].vary_B2c = true;
  mos.mo_ref.opts[0].vary_R0c = {false, true};
  mos.mo_ref.opts[0].vary_Z0s = {false, true};
  mos.mo_ref.opts[0].vary_R0s = {false, false};
  mos.mo_ref.opts[0].vary_Z0c = {false, false};
  mos.mo_ref.opts[0].weight_iota = 30.0;
  mos.mo_ref.opts[0].target_iota = 0.42;
  mos.mo_ref.opts[0].weight_grad_B = 1.0;
  mos.mo_ref.opts[0].weight_B20 = 1.0;
  
  // Set parameters for opt stage 1:
  mos.mo_ref.opts[1].verbose = 0;
  mos.mo_ref.opts[1].fourier_refine = 1;
  mos.mo_ref.opts[1].vary_eta_bar = true;
  mos.mo_ref.opts[1].vary_B2c = true;
  mos.mo_ref.opts[1].vary_R0c = {false, true, true, true};
  mos.mo_ref.opts[1].vary_Z0s = {false, true, true, true};
  mos.mo_ref.opts[1].vary_R0s = {false, false, false, false};
  mos.mo_ref.opts[1].vary_Z0c = {false, false, false, false};
  mos.mo_ref.opts[1].weight_iota = 30.0;
  mos.mo_ref.opts[1].target_iota = 0.42;
  mos.mo_ref.opts[1].weight_grad_B = 1.0;
  mos.mo_ref.opts[1].weight_B20 = 2.0;

  // Run the scan:
  mos.init();
  mos.scan();

  // Check results. Only MPI proc 0 has the final data.
  if (mos.proc0) {
    CHECK(mos.n_scan == 2);

    CHECK(Approx(mos.scan_weight_iota[0]) == 8.0);
    CHECK(Approx(mos.scan_weight_iota[1]) == 10.0);
    CHECK(Approx(mos.scan_weight_B20[0]) == 12.0);
    CHECK(Approx(mos.scan_weight_B20[1]) == 12.0);
    
    CHECK(Approx(mos.scan_R0c(1, 0)) == 0.169146335266187);
    CHECK(Approx(mos.scan_R0c(1, 1)) == 0.17054382023578);
  }
  
  // For each entry in the scan, run a standalone single QSC and
  // confirm all the results match.
  if (mos.proc0) {
    for (int j = 0; j < mos.n_scan; j++) {
      CAPTURE(j);
      
      Qsc q;
      q.nfp = mos.mo_ref.opts[0].q.nfp;
      q.nphi = 25; // Must match value above in this test.
      q.verbose = 0;
      q.eta_bar = mos.scan_eta_bar[j];
      q.sigma0 = mos.scan_sigma0[j];
      q.B2c = mos.scan_B2c[j];
      q.B2s = mos.scan_B2s[j];
      q.order_r_option = "r2.1";
      q.resize_axis_arrays(mos.axis_nmax_plus_1, 0.0);
      for (int k = 0; k < mos.axis_nmax_plus_1; k++) {
	q.R0c[k] = mos.scan_R0c(k, j);
	q.R0s[k] = mos.scan_R0s(k, j);
	q.Z0c[k] = mos.scan_Z0c(k, j);
	q.Z0s[k] = mos.scan_Z0s(k, j);
      }

      q.init();
      q.calculate();

      CHECK(Approx(q.grid_min_R0) == mos.scan_min_R0[j]);
      CHECK(Approx(q.grid_max_curvature) == mos.scan_max_curvature[j]);
      CHECK(Approx(q.iota) == mos.scan_iota[j]);
      CHECK(Approx(q.grid_max_elongation) == mos.scan_max_elongation[j]);
      CHECK(Approx(q.grid_min_L_grad_B) == mos.scan_min_L_grad_B[j]);
      CHECK(Approx(q.grid_min_L_grad_grad_B) == mos.scan_min_L_grad_grad_B[j]);
      CHECK(Approx(q.r_singularity_robust) == mos.scan_r_singularity[j]);
      CHECK(Approx(q.B20_grid_variation) == mos.scan_B20_variation[j]);
      CHECK(Approx(q.B20_residual) == mos.scan_B20_residual[j]);
      CHECK(Approx(q.B20_mean) == mos.scan_B20_mean[j]);
      CHECK(Approx(q.d2_volume_d_psi2) == mos.scan_d2_volume_d_psi2[j]);
      CHECK(Approx(q.DMerc_times_r2) == mos.scan_DMerc_times_r2[j]);
      CHECK(Approx(q.standard_deviation_of_R) == mos.scan_standard_deviation_of_R[j]);
      CHECK(Approx(q.standard_deviation_of_Z) == mos.scan_standard_deviation_of_Z[j]);
      CHECK(q.helicity == mos.scan_helicity[j]);
      CHECK(Approx(q.grid_max_XY2) == mos.scan_max_XY2[j]);
      CHECK(Approx(q.grid_max_Z2) == mos.scan_max_Z2[j]);
      CHECK(Approx(q.grid_max_XY3) == mos.scan_max_XY3[j]);
      CHECK(Approx(q.grid_max_d_XY2_d_varphi) == mos.scan_max_d_XY2_d_varphi[j]);
      CHECK(Approx(q.grid_max_d_Z2_d_varphi) == mos.scan_max_d_Z2_d_varphi[j]);
      CHECK(Approx(q.grid_max_d_XY3_d_varphi) == mos.scan_max_d_XY3_d_varphi[j]);
      CHECK(Approx(q.grid_max_d2_XY2_d_varphi2) == mos.scan_max_d2_XY2_d_varphi2[j]);
      CHECK(Approx(q.grid_max_d2_XY3_d_varphi2) == mos.scan_max_d2_XY3_d_varphi2[j]);
    }
  }
}
