#include <stdexcept>
#include <iostream>
#include "doctest.h"
#include "multiopt.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("Running a standalone opt should yield identical results to a 1-stage multiopt job. [multiopt]") {
  if (single) return;
  int j, k;
  std::string config = "r2 section 5.5";

  for (int vary_axis_option = 0; vary_axis_option < 3; vary_axis_option++) {
    CAPTURE(vary_axis_option);
    for (int vary_scalars_option = 0; vary_scalars_option < 7; vary_scalars_option++) {
      CAPTURE(vary_scalars_option);
      
      Opt opt;
      opt.q = Qsc(config);
      opt.max_iter = 15;

      //opt.q.nphi = 31;
      opt.q.verbose = 0;

      // Turn on all residual terms:
      opt.weight_B20 = 2.0;
      opt.weight_iota = 3.0;
      opt.weight_elongation = 3.5;
      opt.weight_curvature = 3.7;
      opt.weight_R0 = 4.0;
      opt.min_R0 = 0.8;
      opt.weight_d2_volume_d_psi2 = 5.0;
      opt.weight_XY2 = 6.0;
      opt.weight_XY2Prime = 7.0;
      opt.weight_XY2PrimePrime = 7.2;
      opt.weight_Z2 = 6.5;
      opt.weight_Z2Prime = 7.5;
      opt.weight_XY3 = 8.0;
      opt.weight_XY3Prime = 9.0;
      opt.weight_XY3PrimePrime = 9.5;
      opt.weight_grad_B = 10.0;
      opt.weight_grad_grad_B = 11.0;
      opt.weight_r_singularity = 12.0;
      opt.weight_axis_length = 13.0;
      opt.weight_standard_deviation_of_R = 14.0;
      opt.weight_B20_mean = 15.0;
      
      switch (vary_axis_option) {
      case 0:
	// Do not vary the axis
	opt.vary_R0c = {false, false};
	opt.vary_R0s = {false, false};
	opt.vary_Z0c = {false, false};
	opt.vary_Z0s = {false, false};
	break;
      case 1:
	// Do vary the axis, except for the major radius
	opt.vary_R0c = {false, true};
	opt.vary_R0s = {false, true};
	opt.vary_Z0c = {false, true};
	opt.vary_Z0s = {false, true};
	opt.diff_method = DIFF_METHOD_CENTERED;
	break;
      case 2:
	// Vary only selected Fourier modes
	opt.vary_R0c = {false, false};
	opt.vary_R0s = {false, true};
	opt.vary_Z0c = {false, true};
	opt.vary_Z0s = {false, true};
	opt.diff_method = DIFF_METHOD_FORWARD;
	break;
      default:
	CHECK(false); // Should not get here
      }

      switch (vary_scalars_option) {
      case 0:
	opt.vary_eta_bar = true;
	opt.vary_B2c = true;
	opt.vary_B2s = false;
	opt.vary_sigma0 = false;
	break;
      case 1:
	opt.vary_eta_bar = false;
	opt.vary_B2c = false;
	opt.vary_B2s = true;
	opt.vary_sigma0 = true;
	break;
      case 2:
	// Only eta_bar
	opt.vary_eta_bar = true;
	opt.vary_B2c = false;
	opt.vary_B2s = false;
	opt.vary_sigma0 = false;
	break;
      case 3:
	// Only B2c
	opt.vary_eta_bar = false;
	opt.vary_B2c = true;
	opt.vary_B2s = false;
	opt.vary_sigma0 = false;
	break;
      case 4:
	// Only B2s
	opt.vary_eta_bar = false;
	opt.vary_B2c = false;
	opt.vary_B2s = true;
	opt.vary_sigma0 = false;
	break;
      case 5:
	// Only sigma0
	opt.vary_eta_bar = false;
	opt.vary_B2c = false;
	opt.vary_B2s = false;
	opt.vary_sigma0 = true;
	break;
      case 6:
	// Vary everything
	opt.vary_eta_bar = true;
	opt.vary_B2c = true;
	opt.vary_B2s = true;
	opt.vary_sigma0 = true;
	break;
      default:
	CHECK(false); // Should not get here
      }

      // Run the optimization:
      opt.allocate();
      opt.optimize();
      
      // Now set up a 1-stage multiopt with the same parameters:
      MultiOpt mo;
      mo.opts.resize(1);
      mo.opts[0].q = Qsc(config);
      mo.opts[0].q.nphi = opt.q.nphi;
      mo.opts[0].q.verbose = 0;
      mo.opts[0].diff_method = opt.diff_method;
      mo.opts[0].max_iter = opt.max_iter;
      
      mo.opts[0].weight_B20 = opt.weight_B20;
      mo.opts[0].weight_iota = opt.weight_iota;
      mo.opts[0].weight_elongation = opt.weight_elongation;
      mo.opts[0].weight_curvature = opt.weight_curvature;
      mo.opts[0].weight_R0 = opt.weight_R0;
      mo.opts[0].min_R0 = opt.min_R0;
      mo.opts[0].weight_d2_volume_d_psi2 = opt.weight_d2_volume_d_psi2;
      mo.opts[0].weight_XY2 = opt.weight_XY2;
      mo.opts[0].weight_XY2Prime = opt.weight_XY2Prime;
      mo.opts[0].weight_XY2PrimePrime = opt.weight_XY2PrimePrime;
      mo.opts[0].weight_Z2 = opt.weight_Z2;
      mo.opts[0].weight_Z2Prime = opt.weight_Z2Prime;
      mo.opts[0].weight_XY3 = opt.weight_XY3;
      mo.opts[0].weight_XY3Prime = opt.weight_XY3Prime;
      mo.opts[0].weight_XY3PrimePrime = opt.weight_XY3PrimePrime;
      mo.opts[0].weight_grad_B = opt.weight_grad_B;
      mo.opts[0].weight_grad_grad_B = opt.weight_grad_grad_B;
      mo.opts[0].weight_r_singularity = opt.weight_r_singularity;
      mo.opts[0].weight_axis_length = opt.weight_axis_length;
      mo.opts[0].weight_standard_deviation_of_R = opt.weight_standard_deviation_of_R;
      mo.opts[0].weight_B20_mean = opt.weight_B20_mean;
      
      mo.opts[0].vary_R0c = opt.vary_R0c;
      mo.opts[0].vary_R0s = opt.vary_R0s;
      mo.opts[0].vary_Z0c = opt.vary_Z0c;
      mo.opts[0].vary_Z0s = opt.vary_Z0s;
      mo.opts[0].vary_eta_bar = opt.vary_eta_bar;
      mo.opts[0].vary_sigma0 = opt.vary_sigma0;
      mo.opts[0].vary_B2c = opt.vary_B2c;
      mo.opts[0].vary_B2s = opt.vary_B2s;

      // Run the multi-stage optimization:
      mo.optimize();

      REQUIRE(mo.opts[0].n_iter == opt.n_iter);
      CHECK(mo.n_evals == opt.n_evals);
      CHECK(mo.opts[0].n_evals == opt.n_evals);
      
      for (j = 0; j < opt.n_iter; j++) {
	CAPTURE(j);
	// Compare:
	CHECK(Approx(mo.opts[0].iter_eta_bar[j]) == opt.iter_eta_bar[j]);
	CHECK(Approx(mo.opts[0].iter_sigma0[j]) == opt.iter_sigma0[j]);
	CHECK(Approx(mo.opts[0].iter_B2c[j]) == opt.iter_B2c[j]);
	CHECK(Approx(mo.opts[0].iter_B2s[j]) == opt.iter_B2s[j]);
	for (int k = 0; k < 2; k++) {
	  CHECK(Approx(mo.opts[0].iter_R0c(k, j)) == opt.iter_R0c(k, j));
	  CHECK(Approx(mo.opts[0].iter_R0s(k, j)) == opt.iter_R0s(k, j));
	  CHECK(Approx(mo.opts[0].iter_Z0c(k, j)) == opt.iter_Z0c(k, j));
	  CHECK(Approx(mo.opts[0].iter_Z0s(k, j)) == opt.iter_Z0s(k, j));
	}
	CHECK(Approx(mo.opts[0].iter_min_R0[j]) == opt.iter_min_R0[j]);
	CHECK(Approx(mo.opts[0].iter_max_curvature[j]) == opt.iter_max_curvature[j]);
	CHECK(Approx(mo.opts[0].iter_iota[j]) == opt.iter_iota[j]);
	CHECK(Approx(mo.opts[0].iter_max_elongation[j]) == opt.iter_max_elongation[j]);
	CHECK(Approx(mo.opts[0].iter_min_L_grad_B[j]) == opt.iter_min_L_grad_B[j]);
	CHECK(Approx(mo.opts[0].iter_min_L_grad_grad_B[j]) == opt.iter_min_L_grad_grad_B[j]);
	CHECK(Approx(mo.opts[0].iter_r_singularity[j]) == opt.iter_r_singularity[j]);
	CHECK(Approx(mo.opts[0].iter_B20_variation[j]) == opt.iter_B20_variation[j]);
	CHECK(Approx(mo.opts[0].iter_B20_residual[j]) == opt.iter_B20_residual[j]);
	CHECK(Approx(mo.opts[0].iter_d2_volume_d_psi2[j]) == opt.iter_d2_volume_d_psi2[j]);
	CHECK(Approx(mo.opts[0].iter_DMerc_times_r2[j]) == opt.iter_DMerc_times_r2[j]);
	CHECK(Approx(mo.opts[0].iter_standard_deviation_of_R[j]) == opt.iter_standard_deviation_of_R[j]);
	CHECK(Approx(mo.opts[0].iter_standard_deviation_of_Z[j]) == opt.iter_standard_deviation_of_Z[j]);
	CHECK(Approx(mo.opts[0].iter_axis_length[j]) == opt.iter_axis_length[j]);
      
      }
    }
  }
}

TEST_CASE("Check that 2-stage multiopt jobs work for any choice of Fourier refinement. Also check that stage 0 and 1 of the multiopt each match standalone opts with the same parameters. [multiopt]") {
  if (single) return;
  int j, k;
  for (int fourier_refine1 = 0; fourier_refine1 < 4; fourier_refine1++) {
    CAPTURE(fourier_refine1);
    for (int fourier_refine2 = 0; fourier_refine2 < 4; fourier_refine2++) {
      CAPTURE(fourier_refine2);
    
      MultiOpt mo;
      mo.opts.resize(2);

      mo.opts[0].q.verbose = 0;
      mo.opts[0].q.nphi = 31;
      mo.opts[0].q.nfp = 4;
      mo.opts[0].q.R0c = {1.0, 0.17};
      mo.opts[0].q.R0s = {0.0, 0.0};
      mo.opts[0].q.Z0c = {0.0, 0.0};
      mo.opts[0].q.Z0s = {0.0, 0.17};
      mo.opts[0].q.eta_bar = 1.0;
      mo.opts[0].q.order_r_option = "r2.1";
      mo.opts[0].fourier_refine = fourier_refine1;
      mo.opts[0].vary_eta_bar = true;
      mo.opts[0].vary_B2c = true;
      mo.opts[0].vary_R0c = {false, false};
      mo.opts[0].vary_R0s = {false, false};
      mo.opts[0].vary_Z0c = {false, false};
      mo.opts[0].vary_Z0s = {false, true};
      mo.opts[0].weight_grad_B = 1.0;
      mo.opts[0].weight_B20 = 0.1;

      mo.opts[1].q.verbose = 0;
      mo.opts[1].fourier_refine = fourier_refine2;
      mo.opts[1].vary_eta_bar = true;
      mo.opts[1].vary_B2c = true;
      mo.opts[1].weight_grad_B = 1.0;
      mo.opts[1].weight_B20 = 0.3;
      mo.opts[1].weight_grad_grad_B = 1.0e-4;
      mo.opts[1].vary_R0c.resize(2 + fourier_refine1, false);
      mo.opts[1].vary_R0s.resize(2 + fourier_refine1, false);
      mo.opts[1].vary_Z0c.resize(2 + fourier_refine1, false);
      mo.opts[1].vary_Z0s.resize(2 + fourier_refine1, false);
      mo.opts[1].vary_Z0s[1] = true;
      for (j = 0; j < fourier_refine1; j++) {
	mo.opts[1].vary_R0c[j + 2] = true;
	mo.opts[1].vary_Z0s[j + 2] = true;
      }

      // Set up a standalone opt with the same parameters as
      // mo.opts[0].  We do this before running the multiopt because
      // running the multiopt will alter the parameters of the qsc
      // object.
      Opt opt0;
      opt0.q.verbose = 0;
      opt0.q.nphi = mo.opts[0].q.nphi;
      opt0.q.nfp = mo.opts[0].q.nfp;
      opt0.q.eta_bar = mo.opts[0].q.eta_bar;
      opt0.q.sigma0 = mo.opts[0].q.sigma0;
      opt0.q.B2c = mo.opts[0].q.B2c;
      opt0.q.B2s = mo.opts[0].q.B2s;
      opt0.q.order_r_option = mo.opts[0].q.order_r_option;
      opt0.q.R0c = mo.opts[0].q.R0c;
      opt0.q.R0s = mo.opts[0].q.R0s;
      opt0.q.Z0c = mo.opts[0].q.Z0c;
      opt0.q.Z0s = mo.opts[0].q.Z0s;
      
      opt0.diff_method = mo.opts[0].diff_method;
      opt0.max_iter = mo.opts[0].max_iter;
      opt0.fourier_refine = mo.opts[0].fourier_refine;
      opt0.vary_sigma0 = mo.opts[0].vary_sigma0;
      opt0.vary_eta_bar = mo.opts[0].vary_eta_bar;
      opt0.vary_B2c = mo.opts[0].vary_B2c;
      opt0.vary_B2s = mo.opts[0].vary_B2s;
      opt0.vary_R0c = mo.opts[0].vary_R0c;
      opt0.vary_R0s = mo.opts[0].vary_R0s;
      opt0.vary_Z0c = mo.opts[0].vary_Z0c;
      opt0.vary_Z0s = mo.opts[0].vary_Z0s;
      opt0.weight_grad_B = mo.opts[0].weight_grad_B;
      opt0.weight_B20 = mo.opts[0].weight_B20;
      
      mo.optimize();

      int newsize = 2 + fourier_refine1 + fourier_refine2;
      CHECK(mo.opts[1].q.R0c.size() == newsize);
      CHECK(mo.opts[1].q.R0s.size() == newsize);
      CHECK(mo.opts[1].q.Z0c.size() == newsize);
      CHECK(mo.opts[1].q.Z0s.size() == newsize);
      CHECK(mo.opts[1].vary_R0c.size() == newsize);
      CHECK(mo.opts[1].vary_R0s.size() == newsize);
      CHECK(mo.opts[1].vary_Z0c.size() == newsize);
      CHECK(mo.opts[1].vary_Z0s.size() == newsize);

      for (j = 0; j < newsize; j++) {
	CAPTURE(j);
	CHECK(mo.opts[1].vary_R0s[j] == false);
	CHECK(mo.opts[1].vary_Z0c[j] == false);
      }
      for (j = 2; j < newsize; j++) {
	CAPTURE(j);
	CHECK(mo.opts[1].vary_R0c[j]);
	CHECK(mo.opts[1].vary_Z0s[j]);
      }
      CHECK(mo.opts[1].vary_R0c[0] == false);
      CHECK(mo.opts[1].vary_R0c[1] == false);
      CHECK(mo.opts[1].vary_Z0s[0] == false);
      CHECK(mo.opts[1].vary_Z0s[1] == true);

      // The last configuration in stage 0 should match the first configuration in stage 1:
      int index = mo.opts[0].n_iter - 1;
      CHECK(Approx(mo.opts[0].iter_eta_bar[index]) == mo.opts[1].iter_eta_bar[0]);
      CHECK(Approx(mo.opts[0].iter_sigma0[index]) == mo.opts[1].iter_sigma0[0]);
      CHECK(Approx(mo.opts[0].iter_B2c[index]) == mo.opts[1].iter_B2c[0]);
      CHECK(Approx(mo.opts[0].iter_B2s[index]) == mo.opts[1].iter_B2s[0]);
      for (k = 0; k < 2 + fourier_refine1; k++) {
	CHECK(Approx(mo.opts[0].iter_R0c(k, index)) == mo.opts[1].iter_R0c(k, 0));
	CHECK(Approx(mo.opts[0].iter_R0s(k, index)) == mo.opts[1].iter_R0s(k, 0));
	CHECK(Approx(mo.opts[0].iter_Z0c(k, index)) == mo.opts[1].iter_Z0c(k, 0));
	CHECK(Approx(mo.opts[0].iter_Z0s(k, index)) == mo.opts[1].iter_Z0s(k, 0));
      }
      // At the start of stage 1, the new Fourier modes which were not in stage 0 should all be 0:
      for (k = 2 + fourier_refine1; k < newsize; k++) {
	CHECK(Approx(mo.opts[1].iter_R0c(k, 0)) == 0.0);
	CHECK(Approx(mo.opts[1].iter_R0s(k, 0)) == 0.0);
	CHECK(Approx(mo.opts[1].iter_Z0c(k, 0)) == 0.0);
	CHECK(Approx(mo.opts[1].iter_Z0s(k, 0)) == 0.0);
      }

      // Now that the multiopt has run, we know the initial conditions
      // for the latter stage. Set up a standalone opt with the same
      // parameters as mo.opts[1].
      Opt opt1;
      opt1.q.verbose = 0;
      opt1.q.nphi = mo.opts[1].q.nphi;
      opt1.q.nfp = mo.opts[1].q.nfp;
      index = mo.opts[0].n_iter - 1;
      opt1.q.eta_bar = mo.opts[0].iter_eta_bar[index];
      opt1.q.sigma0 = mo.opts[0].iter_sigma0[index];
      opt1.q.B2c = mo.opts[0].iter_B2c[index];
      opt1.q.B2s = mo.opts[0].iter_B2s[index];
      opt1.q.order_r_option = mo.opts[1].q.order_r_option;
      opt1.q.R0c.resize(2 + fourier_refine1);
      opt1.q.R0s.resize(2 + fourier_refine1);
      opt1.q.Z0c.resize(2 + fourier_refine1);
      opt1.q.Z0s.resize(2 + fourier_refine1);
      for (k = 0; k < 2 + fourier_refine1; k++) {
	opt1.q.R0c[k] = mo.opts[0].iter_R0c(k, index);
	opt1.q.R0s[k] = mo.opts[0].iter_R0s(k, index);
	opt1.q.Z0c[k] = mo.opts[0].iter_Z0c(k, index);
	opt1.q.Z0s[k] = mo.opts[0].iter_Z0s(k, index);
      }
      
      opt1.diff_method = mo.opts[1].diff_method;
      opt1.max_iter = mo.opts[1].max_iter;
      opt1.fourier_refine = mo.opts[1].fourier_refine;
      opt1.vary_sigma0 = mo.opts[1].vary_sigma0;
      opt1.vary_eta_bar = mo.opts[1].vary_eta_bar;
      opt1.vary_B2c = mo.opts[1].vary_B2c;
      opt1.vary_B2s = mo.opts[1].vary_B2s;
      opt1.vary_R0c.resize(2 + fourier_refine1, true);
      opt1.vary_R0s.resize(2 + fourier_refine1, false);
      opt1.vary_Z0c.resize(2 + fourier_refine1, false);
      opt1.vary_Z0s.resize(2 + fourier_refine1, true);
      opt1.vary_R0c[0] = false;
      opt1.vary_R0c[1] = false;
      opt1.vary_Z0s[0] = false;
      opt1.weight_grad_B = mo.opts[1].weight_grad_B;
      opt1.weight_B20 = mo.opts[1].weight_B20;
      opt1.weight_grad_grad_B = mo.opts[1].weight_grad_grad_B;
      
      // Run the standalone opts:
      opt0.allocate();
      opt0.optimize();
      opt1.allocate();
      opt1.optimize();

      CHECK(mo.n_evals == opt0.n_evals + opt1.n_evals);
      
      // Compare stage 0 of the multiopt to the standalone opt:
      for (j = 0; j < opt0.n_iter; j++) {
	CAPTURE(j);
	CHECK(Approx(mo.opts[0].iter_eta_bar[j]) == opt0.iter_eta_bar[j]);
	CHECK(Approx(mo.opts[0].iter_sigma0[j]) == opt0.iter_sigma0[j]);
	CHECK(Approx(mo.opts[0].iter_B2c[j]) == opt0.iter_B2c[j]);
	CHECK(Approx(mo.opts[0].iter_B2s[j]) == opt0.iter_B2s[j]);
	for (int k = 0; k < 2 + fourier_refine1; k++) {
	  CHECK(Approx(mo.opts[0].iter_R0c(k, j)) == opt0.iter_R0c(k, j));
	  CHECK(Approx(mo.opts[0].iter_R0s(k, j)) == opt0.iter_R0s(k, j));
	  CHECK(Approx(mo.opts[0].iter_Z0c(k, j)) == opt0.iter_Z0c(k, j));
	  CHECK(Approx(mo.opts[0].iter_Z0s(k, j)) == opt0.iter_Z0s(k, j));
	}
	CHECK(Approx(mo.opts[0].iter_min_R0[j]) == opt0.iter_min_R0[j]);
	CHECK(Approx(mo.opts[0].iter_max_curvature[j]) == opt0.iter_max_curvature[j]);
	CHECK(Approx(mo.opts[0].iter_iota[j]) == opt0.iter_iota[j]);
	CHECK(Approx(mo.opts[0].iter_max_elongation[j]) == opt0.iter_max_elongation[j]);
	CHECK(Approx(mo.opts[0].iter_min_L_grad_B[j]) == opt0.iter_min_L_grad_B[j]);
	CHECK(Approx(mo.opts[0].iter_min_L_grad_grad_B[j]) == opt0.iter_min_L_grad_grad_B[j]);
	CHECK(Approx(mo.opts[0].iter_r_singularity[j]) == opt0.iter_r_singularity[j]);
	CHECK(Approx(mo.opts[0].iter_B20_variation[j]) == opt0.iter_B20_variation[j]);
	CHECK(Approx(mo.opts[0].iter_B20_residual[j]) == opt0.iter_B20_residual[j]);
	CHECK(Approx(mo.opts[0].iter_d2_volume_d_psi2[j]) == opt0.iter_d2_volume_d_psi2[j]);
	CHECK(Approx(mo.opts[0].iter_DMerc_times_r2[j]) == opt0.iter_DMerc_times_r2[j]);
	CHECK(Approx(mo.opts[0].iter_standard_deviation_of_R[j]) == opt0.iter_standard_deviation_of_R[j]);
	CHECK(Approx(mo.opts[0].iter_standard_deviation_of_Z[j]) == opt0.iter_standard_deviation_of_Z[j]);
	CHECK(Approx(mo.opts[0].iter_axis_length[j]) == opt0.iter_axis_length[j]);
      }
      
      // Compare stage 1 of the multiopt to the standalone opt:
      for (j = 0; j < opt1.n_iter; j++) {
	CAPTURE(j);
	CHECK(Approx(mo.opts[1].iter_eta_bar[j]) == opt1.iter_eta_bar[j]);
	CHECK(Approx(mo.opts[1].iter_sigma0[j]) == opt1.iter_sigma0[j]);
	CHECK(Approx(mo.opts[1].iter_B2c[j]) == opt1.iter_B2c[j]);
	CHECK(Approx(mo.opts[1].iter_B2s[j]) == opt1.iter_B2s[j]);
	for (int k = 0; k < newsize; k++) {
	  CHECK(Approx(mo.opts[1].iter_R0c(k, j)) == opt1.iter_R0c(k, j));
	  CHECK(Approx(mo.opts[1].iter_R0s(k, j)) == opt1.iter_R0s(k, j));
	  CHECK(Approx(mo.opts[1].iter_Z0c(k, j)) == opt1.iter_Z0c(k, j));
	  CHECK(Approx(mo.opts[1].iter_Z0s(k, j)) == opt1.iter_Z0s(k, j));
	}
	CHECK(Approx(mo.opts[1].iter_min_R0[j]) == opt1.iter_min_R0[j]);
	CHECK(Approx(mo.opts[1].iter_max_curvature[j]) == opt1.iter_max_curvature[j]);
	CHECK(Approx(mo.opts[1].iter_iota[j]) == opt1.iter_iota[j]);
	CHECK(Approx(mo.opts[1].iter_max_elongation[j]) == opt1.iter_max_elongation[j]);
	CHECK(Approx(mo.opts[1].iter_min_L_grad_B[j]) == opt1.iter_min_L_grad_B[j]);
	CHECK(Approx(mo.opts[1].iter_min_L_grad_grad_B[j]) == opt1.iter_min_L_grad_grad_B[j]);
	CHECK(Approx(mo.opts[1].iter_r_singularity[j]) == opt1.iter_r_singularity[j]);
	CHECK(Approx(mo.opts[1].iter_B20_variation[j]) == opt1.iter_B20_variation[j]);
	CHECK(Approx(mo.opts[1].iter_B20_residual[j]) == opt1.iter_B20_residual[j]);
	CHECK(Approx(mo.opts[1].iter_d2_volume_d_psi2[j]) == opt1.iter_d2_volume_d_psi2[j]);
	CHECK(Approx(mo.opts[1].iter_DMerc_times_r2[j]) == opt1.iter_DMerc_times_r2[j]);
	CHECK(Approx(mo.opts[1].iter_standard_deviation_of_R[j]) == opt1.iter_standard_deviation_of_R[j]);
	CHECK(Approx(mo.opts[1].iter_standard_deviation_of_Z[j]) == opt1.iter_standard_deviation_of_Z[j]);
	CHECK(Approx(mo.opts[1].iter_axis_length[j]) == opt1.iter_axis_length[j]);
      }
      
    }
  }
}

TEST_CASE("Check that using a stage-dependent nphi works. [multiopt]") {
  if (single) return;
  int j, k;
  for (int fourier_refine1 = 1; fourier_refine1 < 2; fourier_refine1++) {
    CAPTURE(fourier_refine1);
    for (int fourier_refine2 = 2; fourier_refine2 < 3; fourier_refine2++) {
      CAPTURE(fourier_refine2);
    
      MultiOpt mo;
      mo.opts.resize(2);

      mo.opts[0].q.verbose = 0;
      mo.opts[0].nphi = {17, 23};
      mo.opts[0].q.nfp = 4;
      mo.opts[0].q.R0c = {1.0, 0.17};
      mo.opts[0].q.R0s = {0.0, 0.0};
      mo.opts[0].q.Z0c = {0.0, 0.0};
      mo.opts[0].q.Z0s = {0.0, 0.17};
      mo.opts[0].q.eta_bar = 1.0;
      mo.opts[0].q.order_r_option = "r2.1";
      mo.opts[0].fourier_refine = fourier_refine1;
      mo.opts[0].vary_eta_bar = true;
      mo.opts[0].vary_B2c = true;
      mo.opts[0].vary_R0c = {false, false};
      mo.opts[0].vary_R0s = {false, false};
      mo.opts[0].vary_Z0c = {false, false};
      mo.opts[0].vary_Z0s = {false, true};
      mo.opts[0].weight_grad_B = 1.0;
      mo.opts[0].weight_B20 = 0.1;

      mo.opts[1].nphi = {27, 31, 37};
      mo.opts[1].q.verbose = 0;
      mo.opts[1].fourier_refine = fourier_refine2;
      mo.opts[1].vary_eta_bar = true;
      mo.opts[1].vary_B2c = true;
      mo.opts[1].weight_grad_B = 1.0;
      mo.opts[1].weight_B20 = 0.3;
      mo.opts[1].weight_grad_grad_B = 1.0e-4;
      mo.opts[1].vary_R0c.resize(2 + fourier_refine1, false);
      mo.opts[1].vary_R0s.resize(2 + fourier_refine1, false);
      mo.opts[1].vary_Z0c.resize(2 + fourier_refine1, false);
      mo.opts[1].vary_Z0s.resize(2 + fourier_refine1, false);
      mo.opts[1].vary_Z0s[1] = true;
      for (j = 0; j < fourier_refine1; j++) {
	mo.opts[1].vary_R0c[j + 2] = true;
	mo.opts[1].vary_Z0s[j + 2] = true;
      }

      // Set up a standalone opt with the same parameters as
      // mo.opts[0].  We do this before running the multiopt because
      // running the multiopt will alter the parameters of the qsc
      // object.
      Opt opt0;
      opt0.q.verbose = 0;
      opt0.nphi = mo.opts[0].nphi;
      opt0.q.nfp = mo.opts[0].q.nfp;
      opt0.q.eta_bar = mo.opts[0].q.eta_bar;
      opt0.q.sigma0 = mo.opts[0].q.sigma0;
      opt0.q.B2c = mo.opts[0].q.B2c;
      opt0.q.B2s = mo.opts[0].q.B2s;
      opt0.q.order_r_option = mo.opts[0].q.order_r_option;
      opt0.q.R0c = mo.opts[0].q.R0c;
      opt0.q.R0s = mo.opts[0].q.R0s;
      opt0.q.Z0c = mo.opts[0].q.Z0c;
      opt0.q.Z0s = mo.opts[0].q.Z0s;
      
      opt0.diff_method = mo.opts[0].diff_method;
      opt0.max_iter = mo.opts[0].max_iter;
      opt0.fourier_refine = mo.opts[0].fourier_refine;
      opt0.vary_sigma0 = mo.opts[0].vary_sigma0;
      opt0.vary_eta_bar = mo.opts[0].vary_eta_bar;
      opt0.vary_B2c = mo.opts[0].vary_B2c;
      opt0.vary_B2s = mo.opts[0].vary_B2s;
      opt0.vary_R0c = mo.opts[0].vary_R0c;
      opt0.vary_R0s = mo.opts[0].vary_R0s;
      opt0.vary_Z0c = mo.opts[0].vary_Z0c;
      opt0.vary_Z0s = mo.opts[0].vary_Z0s;
      opt0.weight_grad_B = mo.opts[0].weight_grad_B;
      opt0.weight_B20 = mo.opts[0].weight_B20;
      
      mo.optimize();

      int newsize = 2 + fourier_refine1 + fourier_refine2;
      CHECK(mo.opts[1].q.R0c.size() == newsize);
      CHECK(mo.opts[1].q.R0s.size() == newsize);
      CHECK(mo.opts[1].q.Z0c.size() == newsize);
      CHECK(mo.opts[1].q.Z0s.size() == newsize);
      CHECK(mo.opts[1].vary_R0c.size() == newsize);
      CHECK(mo.opts[1].vary_R0s.size() == newsize);
      CHECK(mo.opts[1].vary_Z0c.size() == newsize);
      CHECK(mo.opts[1].vary_Z0s.size() == newsize);

      for (j = 0; j < newsize; j++) {
	CAPTURE(j);
	CHECK(mo.opts[1].vary_R0s[j] == false);
	CHECK(mo.opts[1].vary_Z0c[j] == false);
      }
      for (j = 2; j < newsize; j++) {
	CAPTURE(j);
	CHECK(mo.opts[1].vary_R0c[j]);
	CHECK(mo.opts[1].vary_Z0s[j]);
      }
      CHECK(mo.opts[1].vary_R0c[0] == false);
      CHECK(mo.opts[1].vary_R0c[1] == false);
      CHECK(mo.opts[1].vary_Z0s[0] == false);
      CHECK(mo.opts[1].vary_Z0s[1] == true);

      // The last configuration in stage 0 should match the first configuration in stage 1:
      int index = mo.opts[0].n_iter - 1;
      CHECK(Approx(mo.opts[0].iter_eta_bar[index]) == mo.opts[1].iter_eta_bar[0]);
      CHECK(Approx(mo.opts[0].iter_sigma0[index]) == mo.opts[1].iter_sigma0[0]);
      CHECK(Approx(mo.opts[0].iter_B2c[index]) == mo.opts[1].iter_B2c[0]);
      CHECK(Approx(mo.opts[0].iter_B2s[index]) == mo.opts[1].iter_B2s[0]);
      for (k = 0; k < 2 + fourier_refine1; k++) {
	CHECK(Approx(mo.opts[0].iter_R0c(k, index)) == mo.opts[1].iter_R0c(k, 0));
	CHECK(Approx(mo.opts[0].iter_R0s(k, index)) == mo.opts[1].iter_R0s(k, 0));
	CHECK(Approx(mo.opts[0].iter_Z0c(k, index)) == mo.opts[1].iter_Z0c(k, 0));
	CHECK(Approx(mo.opts[0].iter_Z0s(k, index)) == mo.opts[1].iter_Z0s(k, 0));
      }
      // At the start of stage 1, the new Fourier modes which were not in stage 0 should all be 0:
      for (k = 2 + fourier_refine1; k < newsize; k++) {
	CHECK(Approx(mo.opts[1].iter_R0c(k, 0)) == 0.0);
	CHECK(Approx(mo.opts[1].iter_R0s(k, 0)) == 0.0);
	CHECK(Approx(mo.opts[1].iter_Z0c(k, 0)) == 0.0);
	CHECK(Approx(mo.opts[1].iter_Z0s(k, 0)) == 0.0);
      }

      // Now that the multiopt has run, we know the initial conditions
      // for the latter stage. Set up a standalone opt with the same
      // parameters as mo.opts[1].
      Opt opt1;
      opt1.q.verbose = 0;
      opt1.nphi = mo.opts[1].nphi;
      opt1.q.nfp = mo.opts[1].q.nfp;
      index = mo.opts[0].n_iter - 1;
      opt1.q.eta_bar = mo.opts[0].iter_eta_bar[index];
      opt1.q.sigma0 = mo.opts[0].iter_sigma0[index];
      opt1.q.B2c = mo.opts[0].iter_B2c[index];
      opt1.q.B2s = mo.opts[0].iter_B2s[index];
      opt1.q.order_r_option = mo.opts[1].q.order_r_option;
      opt1.q.R0c.resize(2 + fourier_refine1);
      opt1.q.R0s.resize(2 + fourier_refine1);
      opt1.q.Z0c.resize(2 + fourier_refine1);
      opt1.q.Z0s.resize(2 + fourier_refine1);
      for (k = 0; k < 2 + fourier_refine1; k++) {
	opt1.q.R0c[k] = mo.opts[0].iter_R0c(k, index);
	opt1.q.R0s[k] = mo.opts[0].iter_R0s(k, index);
	opt1.q.Z0c[k] = mo.opts[0].iter_Z0c(k, index);
	opt1.q.Z0s[k] = mo.opts[0].iter_Z0s(k, index);
      }
      
      opt1.diff_method = mo.opts[1].diff_method;
      opt1.max_iter = mo.opts[1].max_iter;
      opt1.fourier_refine = mo.opts[1].fourier_refine;
      opt1.vary_sigma0 = mo.opts[1].vary_sigma0;
      opt1.vary_eta_bar = mo.opts[1].vary_eta_bar;
      opt1.vary_B2c = mo.opts[1].vary_B2c;
      opt1.vary_B2s = mo.opts[1].vary_B2s;
      opt1.vary_R0c.resize(2 + fourier_refine1, true);
      opt1.vary_R0s.resize(2 + fourier_refine1, false);
      opt1.vary_Z0c.resize(2 + fourier_refine1, false);
      opt1.vary_Z0s.resize(2 + fourier_refine1, true);
      opt1.vary_R0c[0] = false;
      opt1.vary_R0c[1] = false;
      opt1.vary_Z0s[0] = false;
      opt1.weight_grad_B = mo.opts[1].weight_grad_B;
      opt1.weight_B20 = mo.opts[1].weight_B20;
      opt1.weight_grad_grad_B = mo.opts[1].weight_grad_grad_B;
      
      // Run the standalone opts:
      opt0.allocate();
      opt0.optimize();
      opt1.allocate();
      opt1.optimize();

      CHECK(mo.n_evals == opt0.n_evals + opt1.n_evals);
      
      // Compare stage 0 of the multiopt to the standalone opt:
      for (j = 0; j < opt0.n_iter; j++) {
	CAPTURE(j);
	CHECK(Approx(mo.opts[0].iter_eta_bar[j]) == opt0.iter_eta_bar[j]);
	CHECK(Approx(mo.opts[0].iter_sigma0[j]) == opt0.iter_sigma0[j]);
	CHECK(Approx(mo.opts[0].iter_B2c[j]) == opt0.iter_B2c[j]);
	CHECK(Approx(mo.opts[0].iter_B2s[j]) == opt0.iter_B2s[j]);
	for (int k = 0; k < 2 + fourier_refine1; k++) {
	  CHECK(Approx(mo.opts[0].iter_R0c(k, j)) == opt0.iter_R0c(k, j));
	  CHECK(Approx(mo.opts[0].iter_R0s(k, j)) == opt0.iter_R0s(k, j));
	  CHECK(Approx(mo.opts[0].iter_Z0c(k, j)) == opt0.iter_Z0c(k, j));
	  CHECK(Approx(mo.opts[0].iter_Z0s(k, j)) == opt0.iter_Z0s(k, j));
	}
	CHECK(Approx(mo.opts[0].iter_min_R0[j]) == opt0.iter_min_R0[j]);
	CHECK(Approx(mo.opts[0].iter_max_curvature[j]) == opt0.iter_max_curvature[j]);
	CHECK(Approx(mo.opts[0].iter_iota[j]) == opt0.iter_iota[j]);
	CHECK(Approx(mo.opts[0].iter_max_elongation[j]) == opt0.iter_max_elongation[j]);
	CHECK(Approx(mo.opts[0].iter_min_L_grad_B[j]) == opt0.iter_min_L_grad_B[j]);
	CHECK(Approx(mo.opts[0].iter_min_L_grad_grad_B[j]) == opt0.iter_min_L_grad_grad_B[j]);
	CHECK(Approx(mo.opts[0].iter_r_singularity[j]) == opt0.iter_r_singularity[j]);
	CHECK(Approx(mo.opts[0].iter_B20_variation[j]) == opt0.iter_B20_variation[j]);
	CHECK(Approx(mo.opts[0].iter_B20_residual[j]) == opt0.iter_B20_residual[j]);
	CHECK(Approx(mo.opts[0].iter_d2_volume_d_psi2[j]) == opt0.iter_d2_volume_d_psi2[j]);
	CHECK(Approx(mo.opts[0].iter_DMerc_times_r2[j]) == opt0.iter_DMerc_times_r2[j]);
	CHECK(Approx(mo.opts[0].iter_standard_deviation_of_R[j]) == opt0.iter_standard_deviation_of_R[j]);
	CHECK(Approx(mo.opts[0].iter_standard_deviation_of_Z[j]) == opt0.iter_standard_deviation_of_Z[j]);
	CHECK(Approx(mo.opts[0].iter_axis_length[j]) == opt0.iter_axis_length[j]);
      }
      
      // Compare stage 1 of the multiopt to the standalone opt:
      for (j = 0; j < opt1.n_iter; j++) {
	CAPTURE(j);
	CHECK(Approx(mo.opts[1].iter_eta_bar[j]) == opt1.iter_eta_bar[j]);
	CHECK(Approx(mo.opts[1].iter_sigma0[j]) == opt1.iter_sigma0[j]);
	CHECK(Approx(mo.opts[1].iter_B2c[j]) == opt1.iter_B2c[j]);
	CHECK(Approx(mo.opts[1].iter_B2s[j]) == opt1.iter_B2s[j]);
	for (int k = 0; k < newsize; k++) {
	  CHECK(Approx(mo.opts[1].iter_R0c(k, j)) == opt1.iter_R0c(k, j));
	  CHECK(Approx(mo.opts[1].iter_R0s(k, j)) == opt1.iter_R0s(k, j));
	  CHECK(Approx(mo.opts[1].iter_Z0c(k, j)) == opt1.iter_Z0c(k, j));
	  CHECK(Approx(mo.opts[1].iter_Z0s(k, j)) == opt1.iter_Z0s(k, j));
	}
	CHECK(Approx(mo.opts[1].iter_min_R0[j]) == opt1.iter_min_R0[j]);
	CHECK(Approx(mo.opts[1].iter_max_curvature[j]) == opt1.iter_max_curvature[j]);
	CHECK(Approx(mo.opts[1].iter_iota[j]) == opt1.iter_iota[j]);
	CHECK(Approx(mo.opts[1].iter_max_elongation[j]) == opt1.iter_max_elongation[j]);
	CHECK(Approx(mo.opts[1].iter_min_L_grad_B[j]) == opt1.iter_min_L_grad_B[j]);
	CHECK(Approx(mo.opts[1].iter_min_L_grad_grad_B[j]) == opt1.iter_min_L_grad_grad_B[j]);
	CHECK(Approx(mo.opts[1].iter_r_singularity[j]) == opt1.iter_r_singularity[j]);
	CHECK(Approx(mo.opts[1].iter_B20_variation[j]) == opt1.iter_B20_variation[j]);
	CHECK(Approx(mo.opts[1].iter_B20_residual[j]) == opt1.iter_B20_residual[j]);
	CHECK(Approx(mo.opts[1].iter_d2_volume_d_psi2[j]) == opt1.iter_d2_volume_d_psi2[j]);
	CHECK(Approx(mo.opts[1].iter_DMerc_times_r2[j]) == opt1.iter_DMerc_times_r2[j]);
	CHECK(Approx(mo.opts[1].iter_standard_deviation_of_R[j]) == opt1.iter_standard_deviation_of_R[j]);
	CHECK(Approx(mo.opts[1].iter_standard_deviation_of_Z[j]) == opt1.iter_standard_deviation_of_Z[j]);
	CHECK(Approx(mo.opts[1].iter_axis_length[j]) == opt1.iter_axis_length[j]);
      }
      
    }
  }
}

/*
TEST_CASE("Each stage of a multiopt run should be identical to a standalone opt run. [multiopt]") {
  for (int nopts = 1; nopts < 4; nopts++) {
    CAPTURE(nopts);

    MultiOpt mo;
    mo.opts.resize(nopts);
    mo.verbose = 1;
  }
}
*/
