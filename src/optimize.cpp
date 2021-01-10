#include <fstream>
#include <string>
#include <stdexcept>
#include <cassert>
#include <iomanip>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include "opt.hpp"

int gsl_residual_function(const gsl_vector*, void*, gsl_vector*);
void gsl_callback(const size_t, void*, const gsl_multifit_nlinear_workspace*);

using namespace qsc;

void Opt::optimize() {
  int j;
  if (single) {
    std::cout << "Optimization presently only works with double precision." << std::endl;
    return;
  }
  std::cout << "optimizing..." << std::endl;

  init();
  n_iter = 0;
  q.init();
  
  if (!q.at_least_order_r2) {
    std::cout << "Optimization presently only works with O(r^2)." << std::endl;
    return;
  }

  gsl_vector *gsl_residual = gsl_vector_alloc(n_terms);
  gsl_vector *gsl_state_vector = gsl_vector_alloc(n_parameters);
  gsl_multifit_nlinear_fdf gsl_optimizer;
  gsl_multifit_nlinear_parameters gsl_optimizer_params = gsl_multifit_nlinear_default_parameters();
  gsl_optimizer.f = gsl_residual_function;
  // gsl_optimizer.df = gsl_residual_function_and_Jacobian;
  // gsl_optimizer.fvv = func_fvv;
  gsl_optimizer.df = NULL;
  gsl_optimizer.fvv = NULL;
  gsl_optimizer.n = n_terms;
  gsl_optimizer.p = n_parameters;
  gsl_optimizer.params = (void*)this;
  /*
  // Set finite difference step sizes:
  gsl_optimizer_params.h_df = 1.0e-5;
  gsl_optimizer_params.h_fvv = 1.0e-5;
  */

  // Set initial condition. This would be an invalid cast if SINGLE.
#ifndef SINGLE
  set_state_vector(gsl_state_vector->data);
#endif
  
  switch (algorithm) {
  case GSL_LM:
    gsl_optimizer_params.trs = gsl_multifit_nlinear_trs_lm;
    break;
  case GSL_DOGLEG:
    gsl_optimizer_params.trs = gsl_multifit_nlinear_trs_dogleg;
    break;
  case GSL_DDOGLEG:
    gsl_optimizer_params.trs = gsl_multifit_nlinear_trs_ddogleg;
    break;
  case GSL_SUBSPACE2D:
    gsl_optimizer_params.trs = gsl_multifit_nlinear_trs_subspace2D;
    break;
  default:
    throw std::runtime_error("Error! in optimize_least_squares_gsl.cpp switch! Should not get here!");
  }

  // Set other optimizer parameters
  gsl_optimizer_params.solver = gsl_multifit_nlinear_solver_cholesky;
  // For the above option, there is a trade-off between speed vs
  // robustness when the problem may be rank-deficient. Other options
  // are described in the GSL documentation.
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  const double xtol = 1.0e-8;
  const double gtol = 1.0e-8;
  const double ftol = 1.0e-8;
  gsl_multifit_nlinear_workspace *work = gsl_multifit_nlinear_alloc(T, &gsl_optimizer_params, n_terms, n_parameters);
  gsl_vector * f = gsl_multifit_nlinear_residual(work);
  gsl_vector * x = gsl_multifit_nlinear_position(work);
  int info;

  // GSL runs 1 more iteration than I want it to:
  int max_iter_for_gsl = max_iter - 1;
  
  // Run the optimization
  gsl_multifit_nlinear_init(gsl_state_vector, &gsl_optimizer, work);
  gsl_multifit_nlinear_driver(max_iter_for_gsl, xtol, gtol, ftol,
				gsl_callback, (void*)this, &info, work);

  if (verbose > 0) {
    std::cout << "----- Results from the optimization -----" << std::endl;
    std::cout << "n_iter: " << n_iter << "  niter from GSL: "
	      << gsl_multifit_nlinear_niter(work) << std::endl;
    std::cout << "# of function evals: " << gsl_optimizer.nevalf << std::endl;
  }

  gsl_multifit_nlinear_free(work);
  gsl_vector_free(gsl_residual);
  gsl_vector_free(gsl_state_vector);
  if (verbose > 0) std::cout << "Goodbye from optimize" << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

int gsl_residual_function(const gsl_vector * x, void *params, gsl_vector * f) {
  Opt* opt = (Opt*) params;

  if (opt->verbose > 1) {
    std::cout << "Hello from gsl_residual_function. f stride = " << f->stride << std::endl;
    std::cout << "State vector:";
    for (int j = 0; j < opt->n_parameters; j++)
      std::cout << " " << x->data[j];
    std::cout << std::endl << std::flush;
  }

  /*
  if (n_iter >= max_iter - 1) {
    // Trick GSL into stopping by just returning the same residual as last time,
    // so it thinks the gradient is 0.
    opt->set_residuals(f);
    return GSL_SUCCESS;
  }
  */
  
  // gsl vectors have a 'stride'. Only if the stride is 1 does the layout of a gsl vector correspond to a standard double array.
  // Curran pointed out that the stride for f may not be 1!
  // See https://github.com/PrincetonUniversity/STELLOPT/commit/5820c453283785ffd97e40aec261ca97f76e9071
  assert(x->stride == 1);

  // This unpack_state_vector call would be an invalid cast if SINGLE.
#ifndef SINGLE
  opt->unpack_state_vector(x->data);
#endif
  opt->q.calculate();
  opt->set_residuals(f);

  if (opt->verbose > 1) std::cout << "Goodbye from gsl_residual_function" << std::endl << std::flush;
  return GSL_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////

void Opt::init() {
  int j;
  std::ofstream output_file;

  // Initialize the parameters in the state vector.
  // The order of parameters here must match the order in
  // Opt::set_state_vector() and Opt::unpack_state_vector().
  
  n_parameters = 0;
  if (make_names) state_vector_names.resize(0, "");
  if (vary_eta_bar) {
    if (make_names) state_vector_names.push_back("eta_bar");
    n_parameters++;
  }
  if (vary_sigma0) {
    if (make_names) state_vector_names.push_back("sigma0");
    n_parameters++;
  }
  if (vary_B2c) {
    if (make_names) state_vector_names.push_back("B2c");
    n_parameters++;
  }
  if (vary_B2s) {
    if (make_names) state_vector_names.push_back("B2s");
    n_parameters++;
  }
  for (j = 0; j < vary_R0c.size(); j++) {
    if (vary_R0c[j]) {
      if (make_names) state_vector_names.push_back("R0c[" + std::to_string(j) + "]");
      n_parameters++;
    }
  }
  for (j = 0; j < vary_R0s.size(); j++) {
    if (vary_R0s[j]) {
      if (make_names) state_vector_names.push_back("R0s[" + std::to_string(j) + "]");
      n_parameters++;
    }
  }
  for (j = 0; j < vary_Z0c.size(); j++) {
    if (vary_Z0c[j]) {
      if (make_names) state_vector_names.push_back("Z0c[" + std::to_string(j) + "]");
      n_parameters++;
    }
  }
  for (j = 0; j < vary_Z0s.size(); j++) {
    if (vary_Z0s[j]) {
      if (make_names) state_vector_names.push_back("Z0s[" + std::to_string(j) + "]");
      n_parameters++;
    }
  }

  if (verbose > 0 && make_names) {
    std::cout << "State vector:" << std::endl;
    for (j = 0; j < n_parameters; j++)
      std::cout << "  " << j << ", " << state_vector_names[j] << std::endl;
  }

  // Initialize the terms in the residual.
  // The order of terms here must match the order in
  // Opt::set_residuals().
  n_terms = 0;
  if (make_names) residual_names.resize(0, "");
  if (weight_B20 > 0) {
    n_terms += q.nphi;
    if (make_names)
      for (j = 0; j < q.nphi; j++) residual_names.push_back("B20[" + std::to_string(j) + "]");
  }
  if (weight_iota > 0) {
    n_terms++;
    if (make_names) residual_names.push_back("iota");
  }
  if (weight_R0 > 0) {
    n_terms += q.nphi;
    if (make_names)
      for (j = 0; j < q.nphi; j++) residual_names.push_back("R0[" + std::to_string(j) + "]");
  }
  if (weight_d2_volume_d_psi2 > 0) {
    n_terms++;
    if (make_names) residual_names.push_back("d2_volume_d_psi2");
  }
  if (weight_XY2 > 0) {
    n_terms += q.nphi * 6;
    if (make_names)
      for (j = 0; j < q.nphi * 6; j++) residual_names.push_back("XY2[" + std::to_string(j) + "]");
  }
  if (weight_XY2Prime > 0) {
    n_terms += q.nphi * 6;
    if (make_names)
      for (j = 0; j < q.nphi * 6; j++) residual_names.push_back("XY2Prime[" + std::to_string(j) + "]");
  }
  if (weight_XY3 > 0) {
    n_terms += q.nphi * 8;
    if (make_names)
      for (j = 0; j < q.nphi * 8; j++) residual_names.push_back("XY3[" + std::to_string(j) + "]");
  }
  if (weight_XY3Prime > 0) {
    n_terms += q.nphi * 8;
    if (make_names)
      for (j = 0; j < q.nphi * 8; j++) residual_names.push_back("XY3Prime[" + std::to_string(j) + "]");
  }
  if (weight_grad_grad_B > 0) {
    n_terms += q.nphi * 27;
    if (make_names)
      for (j = 0; j < q.nphi * 27; j++) residual_names.push_back("grad_grad_B[" + std::to_string(j) + "]");
  }
  
  if (make_names) {
    output_file.open("residual_vector");
    if (!output_file.is_open()) {
      throw std::runtime_error("Error! Unable to open output file.");
    }
    assert (n_terms == residual_names.size());
    for (j = 0; j < n_terms; j++)
      output_file << std::setw(4) << j << ", " << residual_names[j] << std::endl;
    output_file.close();
  }

  if (verbose > 0)
    std::cout << "n_parameters: " << n_parameters << "  n_terms: " << n_terms << std::endl;

  if (n_parameters < 1)
    throw std::runtime_error("There must be at least 1 parameter varied.");
  if (n_terms < 1)
    throw std::runtime_error("There must be at least 1 residual term.");
  residuals.resize(n_terms, 0.0);

}

/** Set the optimization state vector from values from the Qsc object.
 */
void Opt::set_state_vector(qscfloat *state_vector) {
  int j, k;
  // The order of parameters here must match the order in Opt::init().
  j = 0;
  if (vary_eta_bar) {
    state_vector[j] = q.eta_bar;
    j++;
  }
  if (vary_sigma0) {
    state_vector[j] = q.sigma0;
    j++;
  }
  if (vary_B2c) {
    state_vector[j] = q.B2c;
    j++;
  }
  if (vary_B2s) {
    state_vector[j] = q.B2s;
    j++;
  }
  for (k = 0; k < vary_R0c.size(); k++) {
    if (vary_R0c[k]) {
      state_vector[j] = q.R0c[k];
      j++;
    }
  }
  for (k = 0; k < vary_R0s.size(); k++) {
    if (vary_R0s[k]) {
      state_vector[j] = q.R0s[k];
      j++;
    }
  }
  for (k = 0; k < vary_Z0c.size(); k++) {
    if (vary_Z0c[k]) {
      state_vector[j] = q.Z0c[k];
      j++;
    }
  }
  for (k = 0; k < vary_Z0s.size(); k++) {
    if (vary_Z0s[k]) {
      state_vector[j] = q.Z0s[k];
      j++;
    }
  }
  assert (j == n_parameters);
}

/** Set parameters of the Qsc object using values from the optimization state vector.
 */
void Opt::unpack_state_vector(qscfloat *state_vector) {
  int j, k;
  // The order of parameters here must match the order in Opt::init().
  j = 0;
  if (vary_eta_bar) {
    q.eta_bar = state_vector[j];
    j++;
  }
  if (vary_sigma0) {
    q.sigma0 = state_vector[j];
    j++;
  }
  if (vary_B2c) {
    q.B2c = state_vector[j];
    j++;
  }
  if (vary_B2s) {
    q.B2s = state_vector[j];
    j++;
  }
  for (k = 0; k < vary_R0c.size(); k++) {
    if (vary_R0c[k]) {
      q.R0c[k] = state_vector[j];
      j++;
    }
  }
  for (k = 0; k < vary_R0s.size(); k++) {
    if (vary_R0s[k]) {
      q.R0s[k] = state_vector[j];
      j++;
    }
  }
  for (k = 0; k < vary_Z0c.size(); k++) {
    if (vary_Z0c[k]) {
      q.Z0c[k] = state_vector[j];
      j++;
    }
  }
  for (k = 0; k < vary_Z0s.size(); k++) {
    if (vary_Z0s[k]) {
      q.Z0s[k] = state_vector[j];
      j++;
    }
  }
  assert (j == n_parameters);
}

/** Set the vector of residuals using values from the Qsc object.
 */
void Opt::set_residuals(gsl_vector* gsl_residual) {
  int j, k;
  j = 0;
  qscfloat term;
  arclength_factor = sqrt(q.d_l_d_phi * (q.d_phi * q.nfp / q.axis_length));

  objective_function = 0.0;
  B20_term = 0.0;
  iota_term = 0.0;
  R0_term = 0.0;
  d2_volume_d_psi2_term = 0.0;
  XY2_term = 0.0;
  XY2Prime_term = 0.0;
  XY3_term = 0.0;
  XY3Prime_term = 0.0;
  grad_grad_B_term = 0.0;

  // The order of terms here must match the order in Opt::init().
  if (weight_B20 > 0) {
    for (k = 0; k < q.nphi; k++) {
      term = weight_B20 * arclength_factor[k] * q.B20_anomaly[k];
      residuals[j] = term;
      B20_term += term * term;
      j++;
    }
  }
  
  if (weight_iota > 0) {
    term = weight_iota - target_iota;
    residuals[j] = term;
    iota_term = term * term;
    j++;
  }

  if (weight_R0 > 0) {
    for (k = 0; k < q.nphi; k++) {
      if (q.R0[k] < min_R0) {
	term = weight_R0 * arclength_factor[k] * (min_R0 - q.R0[k]);
	residuals[j] = term;
	R0_term += term * term;
      } else {
	residuals[j] = 0.0;
      }
      j++;
    }
  }

  if (weight_d2_volume_d_psi2 > 0) {
    if (q.d2_volume_d_psi2 > max_d2_volume_d_psi2) {
      term = weight_d2_volume_d_psi2 * (q.d2_volume_d_psi2 - max_d2_volume_d_psi2);
      residuals[j] = term;
      d2_volume_d_psi2_term = term * term;
    } else {
      residuals[j] = 0.0;
    }
    j++;
  }

  if (weight_XY2 > 0) {
    for (k = 0; k < q.nphi; k++) {
      term = weight_XY2 * arclength_factor[k] * q.X20[k];
      residuals[j] = term;
      XY2_term += term * term;
      j++;
      
      term = weight_XY2 * arclength_factor[k] * q.X2s[k];
      residuals[j] = term;
      XY2_term += term * term;
      j++;
      
      term = weight_XY2 * arclength_factor[k] * q.X2c[k];
      residuals[j] = term;
      XY2_term += term * term;
      j++;
      
      term = weight_XY2 * arclength_factor[k] * q.Y20[k];
      residuals[j] = term;
      XY2_term += term * term;
      j++;
      
      term = weight_XY2 * arclength_factor[k] * q.Y2s[k];
      residuals[j] = term;
      XY2_term += term * term;
      j++;
      
      term = weight_XY2 * arclength_factor[k] * q.Y2c[k];
      residuals[j] = term;
      XY2_term += term * term;
      j++;      
    }
  }

  if (weight_XY2Prime > 0) {
    for (k = 0; k < q.nphi; k++) {
      term = weight_XY2Prime * arclength_factor[k] * q.d_X20_d_varphi[k];
      residuals[j] = term;
      XY2Prime_term += term * term;
      j++;
      
      term = weight_XY2Prime * arclength_factor[k] * q.d_X2s_d_varphi[k];
      residuals[j] = term;
      XY2Prime_term += term * term;
      j++;
      
      term = weight_XY2Prime * arclength_factor[k] * q.d_X2c_d_varphi[k];
      residuals[j] = term;
      XY2Prime_term += term * term;
      j++;
      
      term = weight_XY2Prime * arclength_factor[k] * q.d_Y20_d_varphi[k];
      residuals[j] = term;
      XY2Prime_term += term * term;
      j++;
      
      term = weight_XY2Prime * arclength_factor[k] * q.d_Y2s_d_varphi[k];
      residuals[j] = term;
      XY2Prime_term += term * term;
      j++;
      
      term = weight_XY2Prime * arclength_factor[k] * q.d_Y2c_d_varphi[k];
      residuals[j] = term;
      XY2Prime_term += term * term;
      j++;
    }
  }

  // XY3 and XY3Prime terms should go here.
  if (weight_XY3 > 0) {
    j += q.nphi * 8;
  }
  
  if (weight_XY3Prime > 0) {
    j += q.nphi * 8;
  }
  
  if (weight_grad_grad_B > 0) {
    for (k = 0; k < q.nphi; k++) {
      for (int j2 = 0; j2 < 27; j2++) {
	term = weight_grad_grad_B * arclength_factor[k] * q.grad_grad_B_tensor[k + q.nphi * j2];
	residuals[j] = term;
	grad_grad_B_term += term * term;
	j++;
      }
    }
  }

  // GSL defines the objective function with a factor of 1/2.
  B20_term *= 0.5;
  iota_term *= 0.5;
  R0_term *= 0.5;
  d2_volume_d_psi2_term *= 0.5;
  XY2_term *= 0.5;
  XY2Prime_term *= 0.5;
  XY3_term *= 0.5;
  XY3Prime_term *= 0.5;
  grad_grad_B_term *= 0.5;
    
  objective_function = B20_term + iota_term + R0_term + d2_volume_d_psi2_term
    + XY2_term + XY2Prime_term + XY3_term + XY3Prime_term + grad_grad_B_term;

  assert (j == n_terms);

  // We need to use gsl_vector_set rather than accessing the data of
  // gsl_residual directly because the stride is not 1 when computing
  // finite difference derivatives.
  for (j = 0; j < n_terms; j++)
    gsl_vector_set(gsl_residual, j, residuals[j]);
}

////////////////////////////////////////////////////////////

void gsl_callback(const size_t iter, void *params,
		  const gsl_multifit_nlinear_workspace *w) {
  int j;
  
  Opt* opt = (Opt*) params;
  int n_iter = opt->n_iter; // Shorthand

  gsl_vector *x = gsl_multifit_nlinear_position(w);
  gsl_vector *f = gsl_multifit_nlinear_residual(w);

  if (opt->verbose > 1) {
    std::cout << "gsl_callback called. n_iter=" << opt->n_iter << "  GSL iter=" << iter << std::endl;
    std::cout << "  State vector in callback:";
    for (j = 0; j < opt->n_parameters; j++)
      std::cout << " " << x->data[j];
    std::cout << std::endl << std::flush;
  }

  // A foolproof way to make sure the inputs and outputs are
  // consistent is to just re-run Qsc. This is not as efficient as it
  // could be.
#ifndef SINGLE
  // This next line does not compile if SINGLE.
  opt->unpack_state_vector(x->data);
#endif
  opt->q.calculate();
  opt->set_residuals(f);
  
  opt->iter_objective_function[n_iter] = opt->objective_function;
  opt->iter_B20_term[n_iter] = opt->B20_term;
  opt->iter_iota_term[n_iter] = opt->iota_term;
  opt->iter_R0_term[n_iter] = opt->R0_term;
  opt->iter_d2_volume_d_psi2_term[n_iter] = opt->d2_volume_d_psi2_term;
  opt->iter_XY2_term[n_iter] = opt->XY2_term;
  opt->iter_XY2Prime_term[n_iter] = opt->XY2Prime_term;
  opt->iter_XY3_term[n_iter] = opt->XY3_term;
  opt->iter_XY3Prime_term[n_iter] = opt->XY3Prime_term;
  opt->iter_grad_grad_B_term[n_iter] = opt->grad_grad_B_term;

  opt->iter_eta_bar[n_iter] = opt->q.eta_bar;
  opt->iter_sigma0[n_iter] = opt->q.sigma0;
  opt->iter_B2c[n_iter] = opt->q.B2c;
  opt->iter_B2s[n_iter] = opt->q.B2s;
  opt->iter_min_R0[n_iter] = opt->q.grid_min_R0;
  opt->iter_max_curvature[n_iter] = opt->q.grid_max_curvature;
  opt->iter_iota[n_iter] = opt->q.iota;
  opt->iter_max_elongation[n_iter] = opt->q.grid_max_elongation;
  opt->iter_min_L_grad_B[n_iter] = opt->q.grid_min_L_grad_B;
  opt->iter_min_L_grad_grad_B[n_iter] = opt->q.grid_min_L_grad_grad_B;
  opt->iter_r_singularity[n_iter] = opt->q.r_singularity_robust;
  opt->iter_B20_variation[n_iter] = opt->q.B20_grid_variation;
  opt->iter_B20_residual[n_iter] = opt->q.B20_residual;
  opt->iter_d2_volume_d_psi2[n_iter] = opt->q.d2_volume_d_psi2;
  opt->iter_DMerc_times_r2[n_iter] = opt->q.DMerc_times_r2;
  opt->iter_standard_deviation_of_R[n_iter] = opt->q.standard_deviation_of_R;
  opt->iter_standard_deviation_of_Z[n_iter] = opt->q.standard_deviation_of_Z;
  for (j = 0; j < opt->q.R0c.size(); j++) {
    opt->iter_R0c(j, n_iter) = opt->q.R0c[j];
    opt->iter_R0s(j, n_iter) = opt->q.R0s[j];
    opt->iter_Z0c(j, n_iter) = opt->q.Z0c[j];
    opt->iter_Z0s(j, n_iter) = opt->q.Z0s[j];
  }
    
  opt->n_iter++;
}
