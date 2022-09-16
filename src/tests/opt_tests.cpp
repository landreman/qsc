#include <stdexcept>
#include <iostream>
#include <gsl/gsl_vector.h>
#include "doctest.h"
#include "opt.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("The objective function should equal 1/2 * sum(residuals^2) [opt]") {
  if (single) return;
  for (int j = 0; j < 3; j++) {
    CAPTURE(j);
    
    Opt opt;
    opt.q = Qsc("r2 section 5.5");
    opt.q.verbose = 0;
    opt.q.init();
    opt.q.calculate();

    switch (j) {
    case 0:
      // Turn on all residual terms:
      opt.weight_B20 = 2.0;
      opt.weight_iota = 3.0;
      opt.weight_elongation = 3.5;
      opt.weight_curvature = 3.7;
      opt.weight_R0 = 4.0;
      opt.min_R0 = 0.8;
      opt.weight_d2_volume_d_psi2 = 5.0;
      opt.weight_DMerc_times_r2 = 5.5;
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
      break;
    case 1:
      // Most residual terms off. At least 1 must be on though or qsc will raise an error.
      opt.weight_grad_B = 15.0;
      CHECK(opt.weight_B20 < 0);
      break;
    case 2:
      // Cover the case in which weight_grad_B = 0:
      opt.weight_B20 = 1.0;
      CHECK(opt.weight_grad_B < 0);
      break;
    default:
      CHECK(false); // Should not get here
    }
  
    opt.init_residuals();
    gsl_vector *gsl_residual = gsl_vector_alloc(opt.n_terms);
    opt.set_residuals(gsl_residual);

    Vector temp;
    temp.resize(opt.n_terms, 0.0);
    temp = opt.residuals * opt.residuals;
  
    CHECK(Approx(0.5 * temp.sum()).epsilon(1.0e-13) == opt.objective_function);
  }
}

TEST_CASE("Compute each optimization term a different way and make sure we get the same result. [opt]") {
  if (single) return;

  for (int jconfig = 0; jconfig < 2; jconfig++) {
    CAPTURE(jconfig);
    
    Opt opt;
    if (jconfig == 0) {
      // A config with magnetic well
      opt.q = Qsc("r2 section 5.5");
      opt.min_R0 = 0.98; // so the R0 term is nonzero.
      opt.min_DMerc_times_r2 = 40.0; // so the DMerc term is nonzero.
    } else {
      // A config with magnetic hill
      opt.q = Qsc("r2 section 5.1");
      opt.min_R0 = 0.01; // so R0_term should be 0.
      opt.min_DMerc_times_r2 = 0.0;
    }
    opt.q.verbose = 0;
    opt.q.init();
    opt.q.calculate();

    opt.weight_grad_B = 1.0; // Need at least one residual term.
    opt.init_residuals();
    gsl_vector *gsl_residual = gsl_vector_alloc(opt.n_terms);
    opt.set_residuals(gsl_residual);
  
    qscfloat denominator = opt.q.d_l_d_phi0.sum();
    Vector temp;
    int j;
    qscfloat term;

    temp = opt.q.B20_anomaly * opt.q.B20_anomaly * opt.q.d_l_d_phi0;
    CHECK(Approx(temp.sum() / denominator) == opt.B20_term);

    CHECK(Approx((opt.q.iota - opt.target_iota) * (opt.q.iota - opt.target_iota)) == opt.iota_term);

    temp = opt.q.elongation * opt.q.elongation * opt.q.d_l_d_phi0;
    CHECK(Approx(temp.sum() / denominator) == opt.elongation_term);

    temp = opt.q.curvature * opt.q.curvature * opt.q.d_l_d_phi0;
    CHECK(Approx(temp.sum() / denominator) == opt.curvature_term);

    if (jconfig == 0) {
      CHECK(opt.R0_term > 0);
      temp = (opt.min_R0 - opt.q.R0) * (opt.min_R0 - opt.q.R0) * opt.q.d_l_d_phi0;
      for (j = 0; j < opt.q.nphi; j++) {
	if (opt.q.R0[j] > opt.min_R0) temp[j] = 0;
      }
      CHECK(Approx(temp.sum() / denominator) == opt.R0_term);
    } else {
      CHECK(Approx(opt.R0_term) == 0.0);
    }

    if (jconfig == 0) {
      CHECK(Approx(opt.d2_volume_d_psi2_term) == 0.0);
    } else {
      CHECK(opt.d2_volume_d_psi2_term > 0);
      CHECK(Approx((opt.max_d2_volume_d_psi2 - opt.q.d2_volume_d_psi2) * (opt.max_d2_volume_d_psi2 - opt.q.d2_volume_d_psi2)) == opt.d2_volume_d_psi2_term);
    }

    if (jconfig == 0) {
      CHECK(opt.DMerc_times_r2_term > 0);
      CHECK(Approx((opt.min_DMerc_times_r2 - opt.q.DMerc_times_r2) * (opt.min_DMerc_times_r2 - opt.q.DMerc_times_r2)) == opt.DMerc_times_r2_term);
    } else {
      CHECK(Approx(opt.DMerc_times_r2_term) == 0.0);
    }

    term = 0.0;
    temp = opt.q.X20 * opt.q.X20 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.X2c * opt.q.X2c * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.X2s * opt.q.X2s * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.Y20 * opt.q.Y20 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.Y2c * opt.q.Y2c * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.Y2s * opt.q.Y2s * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    CHECK(Approx(term) == opt.XY2_term);
    
    term = 0.0;
    temp = opt.q.d_X20_d_varphi * opt.q.d_X20_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d_X2c_d_varphi * opt.q.d_X2c_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d_X2s_d_varphi * opt.q.d_X2s_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d_Y20_d_varphi * opt.q.d_Y20_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d_Y2c_d_varphi * opt.q.d_Y2c_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d_Y2s_d_varphi * opt.q.d_Y2s_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    CHECK(Approx(term) == opt.XY2Prime_term);
    
    term = 0.0;
    temp = opt.q.d2_X20_d_varphi2 * opt.q.d2_X20_d_varphi2 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d2_X2c_d_varphi2 * opt.q.d2_X2c_d_varphi2 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d2_X2s_d_varphi2 * opt.q.d2_X2s_d_varphi2 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d2_Y20_d_varphi2 * opt.q.d2_Y20_d_varphi2 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d2_Y2c_d_varphi2 * opt.q.d2_Y2c_d_varphi2 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d2_Y2s_d_varphi2 * opt.q.d2_Y2s_d_varphi2 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    CHECK(Approx(term) == opt.XY2PrimePrime_term);

    term = 0.0;
    temp = opt.q.Z20 * opt.q.Z20 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.Z2c * opt.q.Z2c * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.Z2s * opt.q.Z2s * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    CHECK(Approx(term) == opt.Z2_term);
    
    term = 0.0;
    temp = opt.q.d_Z20_d_varphi * opt.q.d_Z20_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d_Z2c_d_varphi * opt.q.d_Z2c_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d_Z2s_d_varphi * opt.q.d_Z2s_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    CHECK(Approx(term) == opt.Z2Prime_term);

    term = 0.0;
    temp = opt.q.X3c1 * opt.q.X3c1 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.Y3c1 * opt.q.Y3c1 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.Y3s1 * opt.q.Y3s1 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    CHECK(Approx(term) == opt.XY3_term);
    
    term = 0.0;
    temp = opt.q.d_X3c1_d_varphi * opt.q.d_X3c1_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d_Y3c1_d_varphi * opt.q.d_Y3c1_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d_Y3s1_d_varphi * opt.q.d_Y3s1_d_varphi * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    CHECK(Approx(term) == opt.XY3Prime_term);
    
    term = 0.0;
    temp = opt.q.d2_X3c1_d_varphi2 * opt.q.d2_X3c1_d_varphi2 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d2_Y3c1_d_varphi2 * opt.q.d2_Y3c1_d_varphi2 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    temp = opt.q.d2_Y3s1_d_varphi2 * opt.q.d2_Y3s1_d_varphi2 * opt.q.d_l_d_phi0;
    term += temp.sum() / denominator;
    CHECK(Approx(term) == opt.XY3PrimePrime_term);

    term = 0.0;
    for (int j2 = 0; j2 < 3; j2++) {
      for (int j1 = 0; j1 < 3; j1++) {
	for (int jphi = 0; jphi < opt.q.nphi; jphi++) {
	  term += opt.q.grad_B_tensor(jphi, j1, j2) * opt.q.grad_B_tensor(jphi, j1, j2) * opt.q.d_l_d_phi0[jphi];
	}
      }
    }
    CHECK(Approx(term / denominator) == opt.grad_B_term);

    term = 0.0;
    for (int j3 = 0; j3 < 3; j3++) {
      for (int j2 = 0; j2 < 3; j2++) {
	for (int j1 = 0; j1 < 3; j1++) {
	  for (int jphi = 0; jphi < opt.q.nphi; jphi++) {
	    term += opt.q.grad_grad_B_tensor(jphi, j1, j2, j3) * opt.q.grad_grad_B_tensor(jphi, j1, j2, j3) * opt.q.d_l_d_phi0[jphi];
	  }
	}
      }
    }
    CHECK(Approx(term / denominator) == opt.grad_grad_B_term);

    temp = opt.q.d_l_d_phi0 / (opt.q.r_hat_singularity_robust * opt.q.r_hat_singularity_robust);
    CHECK(Approx(temp.sum() / denominator) == opt.r_singularity_term);

    CHECK(Approx(opt.q.axis_length * opt.q.axis_length) == opt.axis_length_term);

    temp = opt.q.d_l_d_phi0 * (opt.q.R0 - opt.q.mean_R) * (opt.q.R0 - opt.q.mean_R);
    CHECK(Approx(temp.sum() / denominator) == opt.standard_deviation_of_R_term);

    temp = opt.q.d_l_d_phi0 * opt.q.B20 * opt.q.B20;
    CHECK(Approx(temp.sum() / denominator) == opt.B20_mean_term);
  }
}

TEST_CASE("Each term in the objective function should be approximately independent of nphi [opt]") {
  if (single) return;
  
  Opt o1, o2;
  o1.q = Qsc("r2 section 5.5");
  o2.q = Qsc("r2 section 5.5");
  // std::cout << "nphi=" << o1.q.nphi << "  R0c=" << o1.q.R0c << " etabar=" << o1.q.eta_bar << " iota=" << o1.q.iota << std::endl;

  o1.q.nphi = 51;
  o2.q.nphi = 101;

  o1.q.init();
  o2.q.init();
  o1.q.calculate();
  o2.q.calculate();

  // Turn on all residual terms:
  o1.weight_B20 = 2.0;
  o1.weight_iota = 3.0;
  o1.weight_elongation = 3.5;
  o1.weight_curvature = 3.7;
  o1.weight_R0 = 4.0;
  o1.min_R0 = 0.8;
  o1.weight_d2_volume_d_psi2 = 5.0;
  o1.weight_DMerc_times_r2 = 5.5;
  o1.weight_XY2 = 6.0;
  o1.weight_XY2Prime = 7.0;
  o1.weight_XY2PrimePrime = 7.2;
  o1.weight_Z2 = 6.5;
  o1.weight_Z2Prime = 7.5;
  o1.weight_XY3 = 8.0;
  o1.weight_XY3Prime = 9.0;
  o1.weight_XY3PrimePrime = 9.5;
  o1.weight_grad_B = 10.0;
  o1.weight_grad_grad_B = 11.0;
  o1.weight_r_singularity = 12.0;
  o1.weight_axis_length = 13.0;
  o1.weight_standard_deviation_of_R = 14.0;
  o1.weight_B20_mean = 15.0;
  
  o2.weight_B20 = 2.0;
  o2.weight_iota = 3.0;
  o2.weight_elongation = 3.5;
  o2.weight_curvature = 3.7;
  o2.weight_R0 = 4.0;
  o2.min_R0 = 0.8;
  o2.weight_d2_volume_d_psi2 = 5.0;
  o2.weight_DMerc_times_r2 = 5.5;
  o2.weight_XY2 = 6.0;
  o2.weight_XY2Prime = 7.0;
  o2.weight_XY2PrimePrime = 7.2;
  o2.weight_Z2 = 6.5;
  o2.weight_Z2Prime = 7.5;
  o2.weight_XY3 = 8.0;
  o2.weight_XY3Prime = 9.0;
  o2.weight_XY3PrimePrime = 9.5;
  o2.weight_grad_B = 10.0;
  o2.weight_grad_grad_B = 11.0;
  o2.weight_r_singularity = 12.0;
  o2.weight_axis_length = 13.0;
  o2.weight_standard_deviation_of_R = 14.0;
  o2.weight_B20_mean = 15.0;

  o1.init_parameters();
  o1.init_residuals();
  o2.init_parameters();
  o2.init_residuals();
  
  gsl_vector *res1 = gsl_vector_alloc(o1.n_terms);
  gsl_vector *res2 = gsl_vector_alloc(o2.n_terms);
  o1.set_residuals(res1);
  o2.set_residuals(res2);

  qscfloat tol = 1.0e-8;
  CHECK(Approx(o1.objective_function).epsilon(tol) == o2.objective_function);
  CHECK(Approx(o1.B20_term).epsilon(tol) == o2.B20_term);
  CHECK(Approx(o1.iota_term).epsilon(tol) == o2.iota_term);
  CHECK(Approx(o1.elongation_term).epsilon(tol) == o2.elongation_term);
  CHECK(Approx(o1.curvature_term).epsilon(tol) == o2.curvature_term);
  CHECK(Approx(o1.R0_term).epsilon(1.0e-4) == o2.R0_term); // This term needs a wider tolerance
  CHECK(Approx(o1.d2_volume_d_psi2_term).epsilon(tol) == o2.d2_volume_d_psi2_term);
  CHECK(Approx(o1.DMerc_times_r2_term).epsilon(tol) == o2.DMerc_times_r2_term);
  CHECK(Approx(o1.XY2_term).epsilon(tol) == o2.XY2_term);
  CHECK(Approx(o1.XY2Prime_term).epsilon(tol) == o2.XY2Prime_term);
  CHECK(Approx(o1.XY2PrimePrime_term).epsilon(tol) == o2.XY2PrimePrime_term);
  CHECK(Approx(o1.Z2_term).epsilon(tol) == o2.Z2_term);
  CHECK(Approx(o1.Z2Prime_term).epsilon(tol) == o2.Z2Prime_term);
  CHECK(Approx(o1.XY3_term).epsilon(tol) == o2.XY3_term);
  CHECK(Approx(o1.XY3Prime_term).epsilon(tol) == o2.XY3Prime_term);
  CHECK(Approx(o1.XY3PrimePrime_term).epsilon(tol) == o2.XY3PrimePrime_term);
  CHECK(Approx(o1.grad_B_term).epsilon(tol) == o2.grad_B_term);
  CHECK(Approx(o1.grad_grad_B_term).epsilon(tol) == o2.grad_grad_B_term);
  CHECK(Approx(o1.r_singularity_term).epsilon(tol) == o2.r_singularity_term);
  CHECK(Approx(o1.axis_length_term).epsilon(tol) == o2.axis_length_term);
  CHECK(Approx(o1.standard_deviation_of_R_term).epsilon(tol) == o2.standard_deviation_of_R_term);
  CHECK(Approx(o1.B20_mean_term).epsilon(tol) == o2.B20_mean_term);

}

TEST_CASE("Running standalone QSC on each configuration in the optimization history should give the corresponding iter_ values. Also confirm that inputs were fixed or varied as requested. [opt]") {
  if (single) return;
  int j, k;

  std::string config = "r2 section 5.5";
  // q0 will be a reference configuration from before the optimization
  Qsc q0(config);
  q0.nphi = 51;
  q0.init();
  q0.calculate();

  for (int vary_axis_option = 0; vary_axis_option < 3; vary_axis_option++) {
    CAPTURE(vary_axis_option);
    for (int vary_scalars_option = 0; vary_scalars_option < 7; vary_scalars_option++) {
      CAPTURE(vary_scalars_option);
      
      Opt opt;
      opt.q = Qsc(config);
      opt.max_iter = 15;

      opt.q.nphi = q0.nphi;
      opt.q.verbose = 0;

      // Turn on all residual terms:
      opt.weight_B20 = 2.0;
      opt.weight_iota = 3.0;
      opt.weight_elongation = 3.5;
      opt.weight_curvature = 3.7;
      opt.weight_R0 = 4.0;
      opt.min_R0 = 0.8;
      opt.weight_d2_volume_d_psi2 = 5.0;
      opt.weight_DMerc_times_r2 = 5.5;
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
      // opt.outfilename = "qsc_out.opt_test.nc";
      // opt.write_netcdf();

      CHECK(opt.n_evals > opt.n_iter);

      // Iteration 0 should exactly match the standalone reference q0:
      CHECK(Approx(opt.iter_eta_bar[0]) == q0.eta_bar);
      CHECK(Approx(opt.iter_B2c[0]) == q0.B2c);
      CHECK(Approx(opt.iter_B2s[0]) == q0.B2s);
      CHECK(Approx(opt.iter_sigma0[0]) == q0.sigma0);
      for (k = 0; k < q0.R0c.size(); k++) {
      	CHECK(Approx(opt.iter_R0c(k, 0)) == q0.R0c[k]);
	CHECK(Approx(opt.iter_R0s(k, 0)) == q0.R0s[k]);
	CHECK(Approx(opt.iter_Z0c(k, 0)) == q0.Z0c[k]);
	CHECK(Approx(opt.iter_Z0s(k, 0)) == q0.Z0s[k]);
      }
      CHECK(Approx(q0.grid_min_R0) == opt.iter_min_R0[0]);
      CHECK(Approx(q0.grid_max_curvature) == opt.iter_max_curvature[0]);
      CHECK(Approx(q0.iota) == opt.iter_iota[0]);
      CHECK(Approx(q0.grid_max_elongation) == opt.iter_max_elongation[0]);
      CHECK(Approx(q0.grid_min_L_grad_B) == opt.iter_min_L_grad_B[0]);
      CHECK(Approx(q0.grid_min_L_grad_grad_B) == opt.iter_min_L_grad_grad_B[0]);
      CHECK(Approx(q0.r_singularity_robust) == opt.iter_r_singularity[0]);
      CHECK(Approx(q0.B20_grid_variation) == opt.iter_B20_variation[0]);
      CHECK(Approx(q0.B20_residual) == opt.iter_B20_residual[0]);
      CHECK(Approx(q0.d2_volume_d_psi2) == opt.iter_d2_volume_d_psi2[0]);
      CHECK(Approx(q0.DMerc_times_r2) == opt.iter_DMerc_times_r2[0]);
      CHECK(Approx(q0.DMerc_times_r2) == opt.iter_DMerc_times_r2[0]);
      CHECK(Approx(q0.standard_deviation_of_R) == opt.iter_standard_deviation_of_R[0]);
      CHECK(Approx(q0.standard_deviation_of_Z) == opt.iter_standard_deviation_of_Z[0]);
      CHECK(Approx(q0.axis_length) == opt.iter_axis_length[0]);

      // Now set up a standalone QSC to check each iteration
      Qsc q;
      q.nphi = opt.q.nphi;
      q.nfp = opt.q.nfp;
      q.I2 = opt.q.I2;
      q.p2 = opt.q.p2;
      q.verbose = 0;
      q.order_r_option = opt.q.order_r_option;
      q.R0c.resize(opt.q.R0c.size(), 0.0);
      q.R0s.resize(opt.q.R0c.size(), 0.0);
      q.Z0c.resize(opt.q.R0c.size(), 0.0);
      q.Z0s.resize(opt.q.R0c.size(), 0.0);
      q.init();
      for (j = 0; j < opt.n_iter; j++) {
	CAPTURE(j);
	// Take inputs from one element of the optimization history and
	// plug them into a standalone QSC:
	q.eta_bar = opt.iter_eta_bar[j];
	q.sigma0 = opt.iter_sigma0[j];
	q.B2c = opt.iter_B2c[j];
	q.B2s = opt.iter_B2s[j];
	for (k = 0; k < opt.q.R0c.size(); k++) {
	  q.R0c[k] = opt.iter_R0c(k, j);
	  q.R0s[k] = opt.iter_R0s(k, j);
	  q.Z0c[k] = opt.iter_Z0c(k, j);
	  q.Z0s[k] = opt.iter_Z0s(k, j);
	}
	// Run standalone QSC:
	q.calculate();
	
	// Compare:
	CHECK(Approx(q.grid_min_R0) == opt.iter_min_R0[j]);
	CHECK(Approx(q.grid_max_curvature) == opt.iter_max_curvature[j]);
	CHECK(Approx(q.iota) == opt.iter_iota[j]);
	CHECK(Approx(q.grid_max_elongation) == opt.iter_max_elongation[j]);
	CHECK(Approx(q.grid_min_L_grad_B) == opt.iter_min_L_grad_B[j]);
	CHECK(Approx(q.grid_min_L_grad_grad_B) == opt.iter_min_L_grad_grad_B[j]);
	CHECK(Approx(q.r_singularity_robust) == opt.iter_r_singularity[j]);
	CHECK(Approx(q.B20_grid_variation) == opt.iter_B20_variation[j]);
	CHECK(Approx(q.B20_residual) == opt.iter_B20_residual[j]);
	CHECK(Approx(q.d2_volume_d_psi2) == opt.iter_d2_volume_d_psi2[j]);
	CHECK(Approx(q.DMerc_times_r2) == opt.iter_DMerc_times_r2[j]);
	CHECK(Approx(q.standard_deviation_of_R) == opt.iter_standard_deviation_of_R[j]);
	CHECK(Approx(q.standard_deviation_of_Z) == opt.iter_standard_deviation_of_Z[j]);
	CHECK(Approx(q.axis_length) == opt.iter_axis_length[j]);
      
	// Check if things were fixed that were supposed to be fixed,
	// and things were varied that were supposed to be varied:
	CHECK(Approx(opt.iter_R0c(0, j)) == q0.R0c[0]);
	CHECK(Approx(opt.iter_R0s(0, j)) == q0.R0s[0]);
	CHECK(Approx(opt.iter_Z0c(0, j)) == q0.Z0c[0]);
	CHECK(Approx(opt.iter_Z0s(0, j)) == q0.Z0s[0]);
	qscfloat eps = 1.0e-13;
	switch (vary_axis_option) {
	case 0:
	  // Do not vary the axis
	  CHECK(Approx(opt.iter_R0c(1, j)) == q0.R0c[1]);
	  CHECK(Approx(opt.iter_R0s(1, j)) == q0.R0s[1]);
	  CHECK(Approx(opt.iter_Z0c(1, j)) == q0.Z0c[1]);
	  CHECK(Approx(opt.iter_Z0s(1, j)) == q0.Z0s[1]);
	  break;
	case 1:
	  // Do vary the axis, except for the major radius
	  if (j > 0) CHECK(Approx(opt.iter_R0c(1, j)).epsilon(eps) != q0.R0c[1]);
	  if (j > 0) CHECK(Approx(opt.iter_R0s(1, j)).epsilon(eps) != q0.R0s[1]);
	  if (j > 0) CHECK(Approx(opt.iter_Z0c(1, j)).epsilon(eps) != q0.Z0c[1]);
	  if (j > 0) CHECK(Approx(opt.iter_Z0s(1, j)).epsilon(eps) != q0.Z0s[1]);
	  break;
	case 2:
	  // Vary only selected Fourier modes
	  CHECK(Approx(opt.iter_R0c(1, j)) == q0.R0c[1]);
	  if (j > 0) CHECK(Approx(opt.iter_R0s(1, j)).epsilon(eps) != q0.R0s[1]);
	  if (j > 0) CHECK(Approx(opt.iter_Z0c(1, j)).epsilon(eps) != q0.Z0c[1]);
	  if (j > 0) CHECK(Approx(opt.iter_Z0s(1, j)).epsilon(eps) != q0.Z0s[1]);
	  break;
	default:
	  CHECK(false); // Should not get here
	}
	
	switch (vary_scalars_option) {
	case 0:
	  if (j > 0) CHECK(Approx(opt.iter_eta_bar[j]) != q0.eta_bar);
	  if (j > 0) CHECK(Approx(opt.iter_B2c[j]) != q0.B2c);
	  CHECK(Approx(opt.iter_B2s[j]) == q0.B2s);
	  CHECK(Approx(opt.iter_sigma0[j]) == q0.sigma0);
	  break;
	case 1:
	  CHECK(Approx(opt.iter_eta_bar[j]) == q0.eta_bar);
	  CHECK(Approx(opt.iter_B2c[j]) == q0.B2c);
	  if (j > 0) CHECK(Approx(opt.iter_B2s[j]) != q0.B2s);
	  if (j > 0) CHECK(Approx(opt.iter_sigma0[j]) != q0.sigma0);
	  break;
	case 2:
	  // Only eta_bar
	  if (j > 0) CHECK(Approx(opt.iter_eta_bar[j]) != q0.eta_bar);
	  CHECK(Approx(opt.iter_B2c[j]) == q0.B2c);
	  CHECK(Approx(opt.iter_B2s[j]) == q0.B2s);
	  CHECK(Approx(opt.iter_sigma0[j]) == q0.sigma0);
	  break;
	case 3:
	  // Only B2c
	  CHECK(Approx(opt.iter_eta_bar[j]) == q0.eta_bar);
	  if (j > 0) CHECK(Approx(opt.iter_B2c[j]) != q0.B2c);
	  CHECK(Approx(opt.iter_B2s[j]) == q0.B2s);
	  CHECK(Approx(opt.iter_sigma0[j]) == q0.sigma0);
	  break;
	case 4:
	  // Only B2s
	  CHECK(Approx(opt.iter_eta_bar[j]) == q0.eta_bar);
	  CHECK(Approx(opt.iter_B2c[j]) == q0.B2c);
	  if (j > 0) CHECK(Approx(opt.iter_B2s[j]) != q0.B2s);
	  CHECK(Approx(opt.iter_sigma0[j]) == q0.sigma0);
	  break;
	case 5:
	  // Only sigma0
	  CHECK(Approx(opt.iter_eta_bar[j]) == q0.eta_bar);
	  CHECK(Approx(opt.iter_B2c[j]) == q0.B2c);
	  CHECK(Approx(opt.iter_B2s[j]) == q0.B2s);
	  if (j > 0) CHECK(Approx(opt.iter_sigma0[j]) != q0.sigma0);
	  break;
	case 6:
	  // Vary everything
	  if (j > 0) CHECK(Approx(opt.iter_eta_bar[j]) != q0.eta_bar);
	  if (j > 0) CHECK(Approx(opt.iter_B2c[j]) != q0.B2c);
	  if (j > 0) CHECK(Approx(opt.iter_B2s[j]) != q0.B2s);
	  if (j > 0) CHECK(Approx(opt.iter_sigma0[j]) != q0.sigma0);
	  break;
	default:
	  CHECK(false); // Should not get here
	}
      }
    }
  }
}

TEST_CASE("Check Opt::unpack_state_vector() and Opt::set_state_vector() [opt]") {
  if (single) return;

  int j, index = 0;
  int n_fourier = 3;
  Vector state_vec1, state_vec2;
  qscfloat arbitrary_val = 3.14;
  Opt opt;
  opt.verbose = 0;
  opt.q.R0c.resize(n_fourier, 0.0);
  opt.q.R0s.resize(n_fourier, 0.0);
  opt.q.Z0c.resize(n_fourier, 0.0);
  opt.q.Z0s.resize(n_fourier, 0.0);
  opt.vary_R0c.resize(n_fourier, false);
  opt.vary_R0s.resize(n_fourier, false);
  opt.vary_Z0c.resize(n_fourier, false);
  opt.vary_Z0s.resize(n_fourier, false);
  opt.weight_grad_B = 1.0; // Need at least 1 residual term.
  // opt.q = Qsc(config);
  std::cout << "max fourier_index:" << (1 << (4 * (n_fourier - 1))) << std::endl;

  for (int vary_eta_bar = 0; vary_eta_bar < 2; vary_eta_bar++) {
    CAPTURE(vary_eta_bar);
    opt.vary_eta_bar = (bool) vary_eta_bar;
    
    for (int vary_sigma0 = 0; vary_sigma0 < 2; vary_sigma0++) {
      CAPTURE(vary_sigma0);
      opt.vary_sigma0 = (bool) vary_sigma0;

      for (int vary_B2c = 0; vary_B2c < 2; vary_B2c++) {
	CAPTURE(vary_B2c);
	opt.vary_B2c = (bool) vary_B2c;

	for (int vary_B2s = 0; vary_B2s < 2; vary_B2s++) {
	  CAPTURE(vary_B2s);
	  opt.vary_B2s = (bool) vary_B2s;

	  for (int vary_R00 = 0; vary_R00 < 2; vary_R00++) {
	    CAPTURE(vary_R00);
	    opt.vary_R0c[0] = (bool) vary_R00;

	    for (int fourier_index = 0; fourier_index < (1 << (4 * (n_fourier - 1))); fourier_index++) {
	      CAPTURE(fourier_index);
	      // std::cout << "fourier_index:" << fourier_index << std::endl;
	      // std::cout << "fourier_index >> (j * 4 + 0):" << (fourier_index >> (j * 4 + 0)) << std::endl;
	      for (j = 0; j < n_fourier - 1; j++) {
		// Pick out the relevant bit of fourier_index:
		opt.vary_R0c[j + 1] = (bool) ((fourier_index >> (j * 4 + 0)) & 1);
		opt.vary_R0s[j + 1] = (bool) ((fourier_index >> (j * 4 + 1)) & 1);
		opt.vary_Z0c[j + 1] = (bool) ((fourier_index >> (j * 4 + 2)) & 1);
		opt.vary_Z0s[j + 1] = (bool) ((fourier_index >> (j * 4 + 3)) & 1);
	      }

	      /*
	      std::cout << "vary_R0c:" << opt.vary_R0c << std::endl;
	      std::cout << "vary_R0s:" << opt.vary_R0s << std::endl;
	      std::cout << "vary_Z0c:" << opt.vary_Z0c << std::endl;
	      std::cout << "vary_Z0s:" << opt.vary_Z0s << std::endl;
	      */
	      index++;
	      if (index == 1) continue; // Need at least 1 parameter to vary
	      opt.init_parameters();
	      opt.init_residuals();
	      state_vec1.resize(opt.n_parameters, 0.0);
	      state_vec2.resize(opt.n_parameters, 0.0);
	      arbitrary_val++;
	      
	      for (j = 0; j < opt.n_parameters; j++) {
		state_vec1[j] = arbitrary_val + 0.1 * j;
	      }
	      opt.unpack_state_vector(&state_vec1[0]);
	      opt.set_state_vector(&state_vec2[0]);
	      for (j = 0; j < opt.n_parameters; j++) {
		CAPTURE(j);
		CHECK(Approx(state_vec2[j]) == state_vec1[j]);
	      }
	    }
	  }
	}
      }
    }
  }
}

TEST_CASE("1d optimization for iota [opt]") {
  if (single) return;

  // Skip the subspace2D algorithm - I'm not sure why it gives an error.
  for (int j_algorithm = 0; j_algorithm < 3; j_algorithm++) {
    CAPTURE(j_algorithm);
    for (int j_var = 0; j_var < 3; j_var++) {
      CAPTURE(j_var);
      Opt opt;
      
      // Set Qsc parameters:
      opt.q.nphi = 31;
      opt.q.nfp = 2;
      opt.q.verbose = 0;
      opt.q.order_r_option = "r2.1";
      opt.q.R0c = {1.0, 0.15};
      opt.q.Z0s = {0.0, 0.15};
      opt.q.R0s = {0.0, 0.0};
      opt.q.Z0c = {0.0, 0.0};
      opt.q.eta_bar = 0.8;
      opt.q.B2c = 0.0;

      // Set optimization parameters:
      opt.verbose = 1;
      switch (j_algorithm) {
      case 0:
	opt.algorithm = GSL_LM;
	break;
      case 1:
	opt.algorithm = GSL_DOGLEG;
	break;
      case 2:
	opt.algorithm = GSL_DDOGLEG;
	break;
      case 3:
	opt.algorithm = GSL_SUBSPACE2D;
	break;
      default:
	CHECK(false); // Should not get here.
      }
      opt.vary_B2c = false;
      switch (j_var) {
      case 0:
	// Vary Z0s[1]
	opt.vary_eta_bar = false;
	opt.vary_R0c = {false, false};
	opt.vary_Z0s = {false, true};
	break;
      case 1:
	// Vary R0c[1]
	opt.vary_eta_bar = false;
	opt.vary_R0c = {false, true};
	opt.vary_Z0s = {false, false};
	opt.diff_method = DIFF_METHOD_CENTERED;
	break;
      case 2:
	// Vary eta_bar
	opt.vary_eta_bar = true;
	opt.vary_R0c = {false, false};
	opt.vary_Z0s = {false, false};
	opt.diff_method = DIFF_METHOD_FORWARD;
	break;
      default:
	CHECK(false); // Should not get here.
      }
      opt.vary_R0s = {false, false};
      opt.vary_Z0c = {false, false};
      
      // Target only iota:
      qscfloat target_iota = -0.5;
      opt.target_iota = target_iota;
      opt.weight_iota = 1.0;
      opt.weight_XY2 = -1.0;
      opt.weight_XY2Prime = -1.0;
      
      opt.allocate();
      opt.optimize();

      CHECK(opt.n_evals > opt.n_iter);
      
      CHECK(Approx(opt.iter_iota[opt.n_iter - 1]) == target_iota);
      switch (j_var) {
      case 0:
	CHECK(Approx(opt.iter_Z0s(1, opt.n_iter - 1)) == 0.0870713164550382);
	break;
      case 1:
	CHECK(Approx(opt.iter_R0c(1, opt.n_iter - 1)) == 0.123975273705787);
	break;
      case 2:
	CHECK(Approx(opt.iter_eta_bar[opt.n_iter - 1]) == 1.03662076014559);
	break;
      default:
	CHECK(false); // Should not get here.
      }
    }
  }
}

TEST_CASE("Try Fourier refinement. Make sure the vary_R0c etc arrays are extended properly. [opt]") {
  if (single) return;

  for (int fourier_refine = 0; fourier_refine < 5; fourier_refine++) {
    CAPTURE(fourier_refine);

    Opt opt;
    // Set up a realistic optimization for QA:
    opt.q.nfp = 2;
    opt.q.nphi = 31;
    opt.q.verbose = 0;
    opt.q.order_r_option = "r2.1";
    opt.q.R0c = {1.0, 0.05};
    opt.q.Z0s = {0.0, -0.05};
    opt.q.R0s = {0.0, 0.0};
    opt.q.Z0c = {0.0, 0.0};
    opt.q.eta_bar = 1.0;
    opt.q.B2c = 0.0;
    opt.vary_eta_bar = true;
    opt.vary_B2c = true;
    opt.vary_R0c = {false, true};
    opt.vary_R0s = {false, false};
    opt.vary_Z0c = {false, false};
    opt.vary_Z0s = {false, true};
    opt.fourier_refine = fourier_refine;
    opt.target_iota = 0.495;
    opt.weight_iota = 1.0e4;
    opt.weight_XY2 = 1.0;
    opt.weight_XY2Prime = 1.0;
    opt.weight_XY3 = 1.0;
    opt.weight_XY3Prime = 1.0;
    opt.weight_B20 = 300.0;
    opt.max_iter = 1000;
    opt.diff_method = DIFF_METHOD_FORWARD;

    opt.allocate();
    opt.optimize();

    CHECK(opt.n_evals > opt.n_iter);
    CHECK(opt.q.nphi == 31);
    REQUIRE(opt.vary_R0c.size() == 2 + fourier_refine);
    REQUIRE(opt.vary_R0s.size() == 2 + fourier_refine);
    REQUIRE(opt.vary_Z0c.size() == 2 + fourier_refine);
    REQUIRE(opt.vary_Z0s.size() == 2 + fourier_refine);
    CHECK_FALSE(opt.vary_R0c[0]);
    CHECK_FALSE(opt.vary_R0s[0]);
    CHECK_FALSE(opt.vary_Z0c[0]);
    CHECK_FALSE(opt.vary_Z0s[0]);
    for (int j = 1; j < 2 + fourier_refine; j++) {
      CAPTURE(j);
      CHECK(opt.vary_R0c[j]);
      CHECK_FALSE(opt.vary_R0s[j]);
      CHECK_FALSE(opt.vary_Z0c[j]);
      CHECK(opt.vary_Z0s[j]);
    }
  }
}

TEST_CASE("Verify that Fourier refinement works gracefully if the max_iter limit is reached [opt]") {
  if (single) return;

  for (int fourier_refine = 0; fourier_refine < 5; fourier_refine++) {
    CAPTURE(fourier_refine);

    Opt opt;
    // Set up a realistic optimization for QA:
    opt.q.nfp = 2;
    opt.q.nphi = 31;
    opt.q.verbose = 0;
    opt.q.order_r_option = "r2.1";
    opt.q.R0c = {1.0, 0.05};
    opt.q.Z0s = {0.0, -0.05};
    opt.q.R0s = {0.0, 0.0};
    opt.q.Z0c = {0.0, 0.0};
    opt.q.eta_bar = 1.0;
    opt.q.B2c = 0.0;
    opt.vary_eta_bar = true;
    opt.vary_B2c = true;
    opt.vary_R0c = {false, true};
    opt.vary_R0s = {false, false};
    opt.vary_Z0c = {false, false};
    opt.vary_Z0s = {false, true};
    opt.fourier_refine = fourier_refine;
    opt.target_iota = 0.495;
    opt.weight_iota = 1.0e4;
    opt.weight_XY2 = 1.0;
    opt.weight_XY2Prime = 1.0;
    opt.weight_XY3 = 1.0;
    opt.weight_XY3Prime = 1.0;
    opt.weight_B20 = 300.0;
    opt.max_iter = 5;
    opt.diff_method = DIFF_METHOD_CENTERED;

    opt.allocate();
    opt.optimize();

    REQUIRE(opt.vary_R0c.size() == 2 + fourier_refine);
    REQUIRE(opt.vary_R0s.size() == 2 + fourier_refine);
    REQUIRE(opt.vary_Z0c.size() == 2 + fourier_refine);
    REQUIRE(opt.vary_Z0s.size() == 2 + fourier_refine);
  }
}

TEST_CASE("Try changing nphi at each stage of Fourier refinement. [opt]") {
  if (single) return;

  int fourier_refine = 2;
  Opt opt;
  // Set up a realistic optimization for QH:
  opt.q.nfp = 4;
  opt.q.nphi = 61; // This value will be over-written by the values in opt.nphi.
  opt.q.verbose = 0;
  opt.q.order_r_option = "r2.1";
  opt.q.R0c = {1.0, 0.17};
  opt.q.Z0s = {0.0, 0.17};
  opt.q.R0s = {0.0, 0.0};
  opt.q.Z0c = {0.0, 0.0};
  opt.q.eta_bar = 1.0;
  opt.q.B2c = 0.0;
  opt.vary_eta_bar = true;
  opt.vary_B2c = true;
  opt.vary_R0c = {false, false};
  opt.vary_R0s = {false, false};
  opt.vary_Z0c = {false, false};
  opt.vary_Z0s = {false, true};
  opt.fourier_refine = fourier_refine;
  opt.nphi = {19, 25, 31};
  opt.weight_grad_B = 1.0;
  opt.weight_B20 = 1.0;

  opt.allocate();
  opt.optimize();

  CHECK(opt.q.nphi == 31);
  REQUIRE(opt.vary_R0c.size() == 2 + fourier_refine);
  REQUIRE(opt.vary_R0s.size() == 2 + fourier_refine);
  REQUIRE(opt.vary_Z0c.size() == 2 + fourier_refine);
  REQUIRE(opt.vary_Z0s.size() == 2 + fourier_refine);
  CHECK_FALSE(opt.vary_R0c[0]);
  CHECK_FALSE(opt.vary_R0s[0]);
  CHECK_FALSE(opt.vary_Z0c[0]);
  CHECK_FALSE(opt.vary_Z0s[0]);
  
  CHECK_FALSE(opt.vary_R0c[1]);
  CHECK_FALSE(opt.vary_R0s[1]);
  CHECK_FALSE(opt.vary_Z0c[1]);
  CHECK(opt.vary_Z0s[1]);
    
  CHECK(opt.vary_R0c[2]);
  CHECK_FALSE(opt.vary_R0s[2]);
  CHECK_FALSE(opt.vary_Z0c[2]);
  CHECK(opt.vary_Z0s[2]);
    
  CHECK(opt.vary_R0c[3]);
  CHECK_FALSE(opt.vary_R0s[3]);
  CHECK_FALSE(opt.vary_Z0c[3]);
  CHECK(opt.vary_Z0s[3]);
}

