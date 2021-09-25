#include <iostream>
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
    CHECK(mos.mo.opts[0].q.nphi == 21);
    CHECK(mos.mo.opts[1].q.nphi == 25);
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

  }
}

TEST_CASE("Run a small MultiOptScan with keep_all false. [mpi] [multiopt_scan]") {
  // This example runs an optimization for QA with nfp = 2.
  
  if (single) return;
  MultiOptScan mos;

  // Set scan parameters:
  mos.verbose = 2;
  mos.keep_all = false;
  mos.min_R0_to_keep = 0.5; // This filter should not be active.
  mos.max_elongation_to_keep = 3.0; // This filter should not be active.
  mos.min_r_singularity_to_keep = 0.25;
  mos.max_B20_variation_to_keep = 0.065;

  mos.params       = {"weight_grad_grad_B", "weight_iota"};
  mos.params_min   = {               0.001,          10.0};
  mos.params_max   = {                 0.1,          30.0};
  mos.params_n     = {                   3,             2};
  mos.params_log   = {                true,          true};
  mos.params_stage = {                   1,            -1};

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
  mos.mo_ref.opts[1].weight_grad_grad_B = 0.01;

  // Run the scan:
  mos.init();
  mos.scan();

  // Check results. Only MPI proc 0 has the final data.
  if (mos.proc0) {
    CHECK(mos.n_scan == 2);

    CHECK(Approx(mos.scan_weight_grad_grad_B[0]) == 0.001);
    CHECK(Approx(mos.scan_weight_grad_grad_B[1]) == 0.001);
    CHECK(Approx(mos.scan_weight_iota[0]) == 10.0);
    CHECK(Approx(mos.scan_weight_iota[1]) == 30.0);
    CHECK(Approx(mos.scan_weight_B20[0]) == 2.0);
    CHECK(Approx(mos.scan_weight_B20[1]) == 2.0);
    
    CHECK(Approx(mos.scan_R0c(1, 0)) == 0.173039969373242);
    CHECK(Approx(mos.scan_R0c(1, 1)) == 0.175276407386024);
  }
}
