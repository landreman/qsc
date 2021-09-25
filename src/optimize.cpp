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

  if (vary_R0c.size() != q.R0c.size()) throw std::runtime_error("Size of vary_R0c is incorrect");
  if (vary_R0s.size() != q.R0s.size()) throw std::runtime_error("Size of vary_R0s is incorrect");
  if (vary_Z0c.size() != q.Z0c.size()) throw std::runtime_error("Size of vary_Z0c is incorrect");
  if (vary_Z0s.size() != q.Z0s.size()) throw std::runtime_error("Size of vary_Z0s is incorrect");
  
  // Add fourier_refine modes to the end of input arrays:
  Vector axis_arr;
  std::valarray<bool> bool_arr;
  int oldsize = q.R0c.size();
  int newsize = oldsize + fourier_refine;
  
  axis_arr = q.R0c;
  q.R0c.resize(newsize, 0.0);
  for (j = 0; j < oldsize; j++) q.R0c[j] = axis_arr[j];
  axis_arr = q.R0s;
  q.R0s.resize(newsize, 0.0);
  for (j = 0; j < oldsize; j++) q.R0s[j] = axis_arr[j];
  axis_arr = q.Z0c;
  q.Z0c.resize(newsize, 0.0);
  for (j = 0; j < oldsize; j++) q.Z0c[j] = axis_arr[j];
  axis_arr = q.Z0s;
  q.Z0s.resize(newsize, 0.0);
  for (j = 0; j < oldsize; j++) q.Z0s[j] = axis_arr[j];
  
  bool_arr = vary_R0c;
  vary_R0c.resize(newsize, false);
  for (j = 0; j < oldsize; j++) vary_R0c[j] = bool_arr[j];
  bool_arr = vary_R0s;
  vary_R0s.resize(newsize, false);
  for (j = 0; j < oldsize; j++) vary_R0s[j] = bool_arr[j];
  bool_arr = vary_Z0c;
  vary_Z0c.resize(newsize, false);
  for (j = 0; j < oldsize; j++) vary_Z0c[j] = bool_arr[j];
  bool_arr = vary_Z0s;
  vary_Z0s.resize(newsize, false);
  for (j = 0; j < oldsize; j++) vary_Z0s[j] = bool_arr[j];

  const gsl_multifit_nlinear_trs * gsl_multifit_nlinear_trs_choice;
  if (verbose > 0) std::cout << "Algorithm: ";
  switch (algorithm) {
  case GSL_LM:
    gsl_multifit_nlinear_trs_choice = gsl_multifit_nlinear_trs_lm;
    if (verbose > 0) std::cout << "Levenberg-Marquardt" << std::endl;
    break;
  case GSL_DOGLEG:
    gsl_multifit_nlinear_trs_choice = gsl_multifit_nlinear_trs_dogleg;
    if (verbose > 0) std::cout << "Dogleg" << std::endl;
    break;
  case GSL_DDOGLEG:
    gsl_multifit_nlinear_trs_choice = gsl_multifit_nlinear_trs_ddogleg;
    if (verbose > 0) std::cout << "Double dogleg" << std::endl;
    break;
  case GSL_SUBSPACE2D:
    gsl_multifit_nlinear_trs_choice = gsl_multifit_nlinear_trs_subspace2D;
    if (verbose > 0) std::cout << "Subspace-2D" << std::endl;
    break;
  default:
    throw std::runtime_error("Error! in optimize_least_squares_gsl.cpp switch! Should not get here!");
  }


  std::cout << "optimizing..." << std::endl;
  
  // init_residuals();
  n_iter = 0;
  for (j_fourier_refine = 0; j_fourier_refine <= fourier_refine; j_fourier_refine++) {
    if (n_iter >= max_iter - 1) {
      if (verbose > 0) std::cout << "Skipping j_fourier_refine = " << j_fourier_refine << " since n_iter is too large." << std::endl;
      break;
    }
    if (j_fourier_refine > 0) {
      vary_R0c[oldsize + j_fourier_refine - 1] = true;
      vary_Z0s[oldsize + j_fourier_refine - 1] = true;
      // For non-stellarator-symmetry, could also modify vary_R0s and vary_Z0c here.
    }
    // If the opt nphi vector is empty, nphi from the original qsc object will be used.
    if (nphi.size() > 0) {
      q.nphi = nphi[j_fourier_refine];
      q.allocate();
      // Changing nphi means the residuals may change size, hence call init_residuals:
      init_residuals();
    } else if (j_fourier_refine == 0) {
      // Make sure init_residuals is always called when j_fourier_refine = 0.
      init_residuals();
    }
    
    if (verbose > 0) {
      std::cout << "`````````````````````````````````````````````````" << std::endl;
      std::cout << "Beginning optimization with j_fourier_refine = " << j_fourier_refine << std::endl;
      std::cout << "vary_R0c: ";
      std::cout << vary_R0c << std::endl;
      std::cout << "vary_R0s: ";
      std::cout << vary_R0s << std::endl;
      std::cout << "vary_Z0c: ";
      std::cout << vary_Z0c << std::endl;
      std::cout << "vary_Z0s: ";
      std::cout << vary_Z0s << std::endl;
      std::cout << "nphi: " << q.nphi << std::endl;
    }
    init_parameters();
    q.init();
  
    if (!q.order_r2p1) {
      std::cout << "Optimization presently only works with order_r_option = r2.1." << std::endl;
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

    // Set other optimizer parameters
    gsl_optimizer_params.trs = gsl_multifit_nlinear_trs_choice;
    gsl_optimizer_params.solver = gsl_multifit_nlinear_solver_cholesky;
    // For the above option, there is a trade-off between speed vs
    // robustness when the problem may be rank-deficient. Other options
    // are described in the GSL documentation.
    switch (diff_method) {
    case DIFF_METHOD_FORWARD:
      gsl_optimizer_params.fdtype = GSL_MULTIFIT_NLINEAR_FWDIFF;
      break;
    case DIFF_METHOD_CENTERED:
      gsl_optimizer_params.fdtype = GSL_MULTIFIT_NLINEAR_CTRDIFF;
      break;
    default:
      throw std::runtime_error("Unrecognized diff_method.");
    }
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    const double xtol = 1.0e-8;
    const double gtol = 1.0e-8;
    const double ftol = 1.0e-8;
    gsl_multifit_nlinear_workspace *work = gsl_multifit_nlinear_alloc(T, &gsl_optimizer_params, n_terms, n_parameters);
    gsl_vector * f = gsl_multifit_nlinear_residual(work);
    gsl_vector * x = gsl_multifit_nlinear_position(work);
    int info;
    
    // GSL runs 1 more iteration than I want it to:
    int max_iter_for_gsl = max_iter - 1 - n_iter;
    
    // Run the optimization
    gsl_multifit_nlinear_init(gsl_state_vector, &gsl_optimizer, work);
    gsl_multifit_nlinear_driver(max_iter_for_gsl, xtol, gtol, ftol,
				gsl_callback, (void*)this, &info, work);
    
    if (verbose > 0) {
      std::cout << "----- Results from the optimization -----" << std::endl;
      std::cout << "n_iter: " << n_iter << "  niter from GSL: "
		<< gsl_multifit_nlinear_niter(work) << std::endl;
      std::cout << "# of function evals: " << gsl_optimizer.nevalf << std::endl;
      std::cout << "Final configuration:" << std::endl;
      std::cout << "  eta_bar: " << q.eta_bar << "  sigma0: " << q.sigma0 << std::endl;
      std::cout << "  B2c: " << q.B2c << "  B2s: " << q.B2s << std::endl;
      std::cout << "  R0c: " << q.R0c << std::endl;
      std::cout << "  R0s: " << q.R0s << std::endl;
      std::cout << "  Z0c: " << q.Z0c << std::endl;
      std::cout << "  Z0s: " << q.Z0s << std::endl;
      std::cout << "  iota: " << q.iota << "  r_singularity: " << q.r_singularity_robust << std::endl;
      std::cout << "  L grad B: " << q.grid_min_L_grad_B << "  L grad grad B: " << q.grid_min_L_grad_grad_B << std::endl;
      std::cout << "  B20 variation: " << q.B20_grid_variation << "  min(R0): " << q.grid_min_R0 << std::endl;
      std::cout << "  d2 volume / d psi2: " << q.d2_volume_d_psi2 << std::endl;
      std::cout << "  axis_length: " << q.axis_length << "  stddev(R): " << q.standard_deviation_of_R << std::endl;
    }
    
    gsl_multifit_nlinear_free(work);
    gsl_vector_free(gsl_residual);
    gsl_vector_free(gsl_state_vector);
  }
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

void Opt::init_parameters() {
  int j;
 
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
    std::cout << "State vector names:" << std::endl;
    for (j = 0; j < n_parameters; j++)
      std::cout << "  " << j << ", " << state_vector_names[j] << std::endl;
  }
  
  if (verbose > 0)
    std::cout << "n_parameters: " << n_parameters << std::endl;

  if (n_parameters < 1)
    throw std::runtime_error("There must be at least 1 parameter varied.");
}

//////////////////////////////////////////////////////////////////////////////

void Opt::init_residuals() {
  int j;
  std::ofstream output_file;

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
  if (weight_elongation > 0) {
    n_terms += q.nphi;
    if (make_names)
      for (j = 0; j < q.nphi; j++) residual_names.push_back("elongation[" + std::to_string(j) + "]");
  }
  if (weight_curvature > 0) {
    n_terms += q.nphi;
    if (make_names)
      for (j = 0; j < q.nphi; j++) residual_names.push_back("curvature[" + std::to_string(j) + "]");
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
  if (weight_Z2 > 0) {
    n_terms += q.nphi * 3;
    if (make_names)
      for (j = 0; j < q.nphi * 3; j++) residual_names.push_back("Z2[" + std::to_string(j) + "]");
  }
  if (weight_Z2Prime > 0) {
    n_terms += q.nphi * 3;
    if (make_names)
      for (j = 0; j < q.nphi * 3; j++) residual_names.push_back("Z2Prime[" + std::to_string(j) + "]");
  }
  if (weight_XY3 > 0) {
    n_terms += q.nphi * 3;
    if (make_names)
      for (j = 0; j < q.nphi * 3; j++) residual_names.push_back("XY3[" + std::to_string(j) + "]");
  }
  if (weight_XY3Prime > 0) {
    n_terms += q.nphi * 3;
    if (make_names)
      for (j = 0; j < q.nphi * 3; j++) residual_names.push_back("XY3Prime[" + std::to_string(j) + "]");
  }
  if (weight_grad_B > 0) {
    n_terms += q.nphi * 9;
    if (make_names)
      for (j = 0; j < q.nphi * 9; j++) residual_names.push_back("grad_B[" + std::to_string(j) + "]");
  }
  if (weight_grad_grad_B > 0) {
    n_terms += q.nphi * 27;
    if (make_names)
      for (j = 0; j < q.nphi * 27; j++) residual_names.push_back("grad_grad_B[" + std::to_string(j) + "]");
  }
  if (weight_r_singularity > 0) {
    n_terms += q.nphi;
    if (make_names)
      for (j = 0; j < q.nphi; j++) residual_names.push_back("r_singularity[" + std::to_string(j) + "]");
  }
  if (weight_axis_length > 0) {
    n_terms += 1;
    if (make_names) residual_names.push_back("axis_length");
  }
  if (weight_standard_deviation_of_R > 0) {
    n_terms += q.nphi;
    if (make_names)
      for (j = 0; j < q.nphi; j++) residual_names.push_back("standard_deviation_of_R[" + std::to_string(j) + "]");
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
    std::cout << "n_terms: " << n_terms << std::endl;

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
  elongation_term = 0.0;
  curvature_term = 0.0;
  R0_term = 0.0;
  d2_volume_d_psi2_term = 0.0;
  XY2_term = 0.0;
  XY2Prime_term = 0.0;
  Z2_term = 0.0;
  Z2Prime_term = 0.0;
  XY3_term = 0.0;
  XY3Prime_term = 0.0;
  grad_B_term = 0.0;
  grad_grad_B_term = 0.0;
  r_singularity_term = 0.0;
  axis_length_term = 0.0;
  standard_deviation_of_R_term = 0.0;

  // The order of terms here must match the order in Opt::init_residuals().
  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * q.B20_anomaly[k];
    B20_term += term * term;
    if (weight_B20 > 0) residuals[j++] = weight_B20 * term;
  }
  
  term = q.iota - target_iota;
  iota_term = term * term;
  if (weight_iota > 0) residuals[j++] = weight_iota * term;

  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * q.elongation[k];
    elongation_term += term * term;
    if (weight_elongation > 0) residuals[j++] = weight_elongation * term;
  }
  
  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * q.curvature[k];
    curvature_term += term * term;
    if (weight_curvature > 0) residuals[j++] = weight_curvature * term;
  }
  
  for (k = 0; k < q.nphi; k++) {
    if (q.R0[k] < min_R0) {
      term = arclength_factor[k] * (min_R0 - q.R0[k]);
    } else {
      term = 0.0;
    }
    R0_term += term * term;
    if (weight_R0 > 0) residuals[j++] = weight_R0 * term;
  }

  if (q.d2_volume_d_psi2 > max_d2_volume_d_psi2) {
    term = q.d2_volume_d_psi2 - max_d2_volume_d_psi2;
  } else {
    term = 0.0;
  }
  if (weight_d2_volume_d_psi2 > 0) residuals[j++] = weight_d2_volume_d_psi2 * term;
  d2_volume_d_psi2_term = term * term;

  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * q.X20[k];
    XY2_term += term * term;
    if (weight_XY2 > 0) residuals[j++] = weight_XY2 * term;
    
    term = arclength_factor[k] * q.X2s[k];
    XY2_term += term * term;
    if (weight_XY2 > 0) residuals[j++] = weight_XY2 * term;
    
    term = arclength_factor[k] * q.X2c[k];
    XY2_term += term * term;
    if (weight_XY2 > 0) residuals[j++] = weight_XY2 * term;
    
    term = arclength_factor[k] * q.Y20[k];
    XY2_term += term * term;
    if (weight_XY2 > 0) residuals[j++] = weight_XY2 * term;
    
    term = arclength_factor[k] * q.Y2s[k];
    XY2_term += term * term;
    if (weight_XY2 > 0) residuals[j++] = weight_XY2 * term;
    
    term = arclength_factor[k] * q.Y2c[k];
    XY2_term += term * term;
    if (weight_XY2 > 0) residuals[j++] = weight_XY2 * term;
  }

  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * q.d_X20_d_varphi[k];
    XY2Prime_term += term * term;
    if (weight_XY2Prime > 0) residuals[j++] = weight_XY2Prime * term;
    
    term = arclength_factor[k] * q.d_X2s_d_varphi[k];
    XY2Prime_term += term * term;
    if (weight_XY2Prime > 0) residuals[j++] = weight_XY2Prime * term;
    
    term = arclength_factor[k] * q.d_X2c_d_varphi[k];
    XY2Prime_term += term * term;
    if (weight_XY2Prime > 0) residuals[j++] = weight_XY2Prime * term;
    
    term = arclength_factor[k] * q.d_Y20_d_varphi[k];
    XY2Prime_term += term * term;
    if (weight_XY2Prime > 0) residuals[j++] = weight_XY2Prime * term;
    
    term = arclength_factor[k] * q.d_Y2s_d_varphi[k];
    XY2Prime_term += term * term;
    if (weight_XY2Prime > 0) residuals[j++] = weight_XY2Prime * term;
    
    term = arclength_factor[k] * q.d_Y2c_d_varphi[k];
    XY2Prime_term += term * term;
    if (weight_XY2Prime > 0) residuals[j++] = weight_XY2Prime * term;
  }

  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * q.Z20[k];
    Z2_term += term * term;
    if (weight_Z2 > 0) residuals[j++] = weight_Z2 * term;
    
    term = arclength_factor[k] * q.Z2s[k];
    Z2_term += term * term;
    if (weight_Z2 > 0) residuals[j++] = weight_Z2 * term;
    
    term = arclength_factor[k] * q.Z2c[k];
    Z2_term += term * term;
    if (weight_Z2 > 0) residuals[j++] = weight_Z2 * term;
  }

  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * q.d_Z20_d_varphi[k];
    Z2Prime_term += term * term;
    if (weight_Z2Prime > 0) residuals[j++] = weight_Z2Prime * term;
      
    term = arclength_factor[k] * q.d_Z2s_d_varphi[k];
    Z2Prime_term += term * term;
    if (weight_Z2Prime > 0) residuals[j++] = weight_Z2Prime * term;
      
    term = arclength_factor[k] * q.d_Z2c_d_varphi[k];
    Z2Prime_term += term * term;
    if (weight_Z2Prime > 0) residuals[j++] = weight_Z2Prime * term;
  }

  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * q.X3c1[k];
    XY3_term += term * term;
    if (weight_XY3 > 0) residuals[j++] = weight_XY3 * term;
      
    term = arclength_factor[k] * q.Y3c1[k];
    XY3_term += term * term;
    if (weight_XY3 > 0) residuals[j++] = weight_XY3 * term;
      
    term = arclength_factor[k] * q.Y3s1[k];
    XY3_term += term * term;
    if (weight_XY3 > 0) residuals[j++] = weight_XY3 * term;
  }
  
  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * q.d_X3c1_d_varphi[k];
    XY3Prime_term += term * term;
    if (weight_XY3Prime > 0) residuals[j++] = weight_XY3Prime * term;
      
    term = arclength_factor[k] * q.d_Y3c1_d_varphi[k];
    XY3Prime_term += term * term;
    if (weight_XY3Prime > 0) residuals[j++] = weight_XY3Prime * term;
      
    term = arclength_factor[k] * q.d_Y3s1_d_varphi[k];
    XY3Prime_term += term * term;
    if (weight_XY3Prime > 0) residuals[j++] = weight_XY3Prime * term;
  }
  
  for (k = 0; k < q.nphi; k++) {
    for (int j2 = 0; j2 < 9; j2++) {
      term = arclength_factor[k] * q.grad_B_tensor[k + q.nphi * j2];
      grad_B_term += term * term;
      if (weight_grad_B > 0) residuals[j++] = weight_grad_B * term;
    }
  }

  for (k = 0; k < q.nphi; k++) {
    for (int j2 = 0; j2 < 27; j2++) {
      term = arclength_factor[k] * q.grad_grad_B_tensor[k + q.nphi * j2];
      grad_grad_B_term += term * term;
      if (weight_grad_grad_B > 0) residuals[j++] = weight_grad_grad_B * term;
    }
  }

  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] / q.r_hat_singularity_robust[k];
    r_singularity_term += term * term;
    if (weight_r_singularity > 0) residuals[j++] = weight_r_singularity * term;
  }

  term = q.axis_length;
  axis_length_term = term * term;
  if (weight_axis_length > 0) residuals[j++] = weight_axis_length * term;

  for (k = 0; k < q.nphi; k++) {
    term = arclength_factor[k] * (q.R0[k] - q.mean_R);
    standard_deviation_of_R_term += term * term;
    if (weight_standard_deviation_of_R > 0) residuals[j++] = weight_standard_deviation_of_R * term;
  }

  Vector weights = {weight_B20, weight_iota, weight_elongation, weight_curvature, weight_R0, weight_d2_volume_d_psi2,
    weight_XY2, weight_XY2Prime, weight_Z2, weight_Z2Prime, weight_XY3, weight_XY3Prime,
    weight_grad_B, weight_grad_grad_B, weight_r_singularity,
    weight_axis_length, weight_standard_deviation_of_R};

  Vector terms = {B20_term, iota_term, elongation_term, curvature_term,
    R0_term, d2_volume_d_psi2_term, XY2_term, XY2Prime_term, Z2_term, Z2Prime_term,
    XY3_term, XY3Prime_term, grad_B_term, grad_grad_B_term, r_singularity_term,
    axis_length_term, standard_deviation_of_R_term};

  assert (weights.size() == terms.size());
  for (k = 0; k < terms.size(); k++) {
    // GSL defines the objective function with a factor of 1/2.
    if (weights[k] > 0) objective_function += 0.5 * weights[k] * weights[k] * terms[k];
  }

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
  opt->iter_elongation_term[n_iter] = opt->elongation_term;
  opt->iter_curvature_term[n_iter] = opt->curvature_term;
  opt->iter_R0_term[n_iter] = opt->R0_term;
  opt->iter_d2_volume_d_psi2_term[n_iter] = opt->d2_volume_d_psi2_term;
  opt->iter_XY2_term[n_iter] = opt->XY2_term;
  opt->iter_XY2Prime_term[n_iter] = opt->XY2Prime_term;
  opt->iter_Z2_term[n_iter] = opt->Z2_term;
  opt->iter_Z2Prime_term[n_iter] = opt->Z2Prime_term;
  opt->iter_XY3_term[n_iter] = opt->XY3_term;
  opt->iter_XY3Prime_term[n_iter] = opt->XY3Prime_term;
  opt->iter_grad_B_term[n_iter] = opt->grad_B_term;
  opt->iter_grad_grad_B_term[n_iter] = opt->grad_grad_B_term;
  opt->iter_r_singularity_term[n_iter] = opt->r_singularity_term;
  opt->iter_axis_length_term[n_iter] = opt->axis_length_term;
  opt->iter_standard_deviation_of_R_term[n_iter] = opt->standard_deviation_of_R_term;

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
  opt->iter_axis_length[n_iter] = opt->q.axis_length;
  for (j = 0; j < opt->q.R0c.size(); j++) {
    opt->iter_R0c(j, n_iter) = opt->q.R0c[j];
    opt->iter_R0s(j, n_iter) = opt->q.R0s[j];
    opt->iter_Z0c(j, n_iter) = opt->q.Z0c[j];
    opt->iter_Z0s(j, n_iter) = opt->q.Z0s[j];
  }

  opt->iter_fourier_refine_step[n_iter] = opt->j_fourier_refine;
  opt->n_iter++;
}
