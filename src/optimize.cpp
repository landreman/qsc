#include <fstream>
#include <string>
#include <stdexcept>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include "opt.hpp"

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
  gsl_vector *gsl_residual = gsl_vector_alloc(n_terms);
  gsl_vector *gsl_state_vector = gsl_vector_alloc(n_parameters);
  gsl_multifit_nlinear_fdf gsl_optimizer;
  gsl_multifit_nlinear_parameters gsl_optimizer_params = gsl_multifit_nlinear_default_parameters();
}

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
    n_terms += q.nphi * 6;
    if (make_names)
      for (j = 0; j < q.nphi * 6; j++) residual_names.push_back("XY3[" + std::to_string(j) + "]");
  }
  if (weight_XY3Prime > 0) {
    n_terms += q.nphi * 6;
    if (make_names)
      for (j = 0; j < q.nphi * 6; j++) residual_names.push_back("XY3Prime[" + std::to_string(j) + "]");
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
void Opt::set_residuals(qscfloat* residuals) {
  int j, k;
  j = 0;
  qscfloat term;
  arclength_factor = sqrt(q.d_l_d_phi) * (q.d_phi * q.nfp / q.axis_length);

  qscfloat objective_function = 0.0;
  qscfloat B20_term = 0.0, iota_term = 0.0;
  qscfloat R0_term = 0.0, d2_volume_d_psi2_term = 0.0;
  qscfloat XY2_term = 0.0, XY2Prime_term = 0.0;
  qscfloat XY3_term = 0.0, XY3Prime_term = 0.0;
  qscfloat grad_grad_B_term = 0.0;

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

  if (weight_grad_grad_B > 0) {
    for (k = 0; k < q.nphi * 27; k++) {
      term = weight_grad_grad_B * arclength_factor[k] * q.grad_grad_B_tensor[k];
      residuals[j] = term;
      grad_grad_B_term += term * term;
      j++;
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

  iter_objective_function[n_iter] = objective_function;
  iter_B20_term[n_iter] = B20_term;
  iter_iota_term[n_iter] = iota_term;
  iter_R0_term[n_iter] = R0_term;
  iter_d2_volume_d_psi2_term[n_iter] = d2_volume_d_psi2_term;
  iter_XY2_term[n_iter] = XY2_term;
  iter_XY2Prime_term[n_iter] = XY2Prime_term;
  iter_XY3_term[n_iter] = XY3_term;
  iter_XY3Prime_term[n_iter] = XY3Prime_term;
  iter_grad_grad_B_term[n_iter] = grad_grad_B_term;
  
  assert (j == n_terms);
}
