#include <iostream>
#include <valarray>
#include "qsc.hpp"

using namespace qsc;

/** Residual for the sigma equation.
 */
void Qsc::sigma_eq_residual(Vector& state, Vector& residual, void* user_data) {
  // Get a pointer to the relevant Qsc object:
  Qsc* q = (Qsc*)user_data;

  q->sigma = state;
  q->sigma[0] = q->sigma0;
  q->iota = state[0];

  // residual = d_d_varphi * sigma;
  matrix_vector_product(q->d_d_varphi, q->sigma, q->residual);

  q->residual += (q->iota + q->helicity * q->nfp)
    * (q->etabar_squared_over_curvature_squared * q->etabar_squared_over_curvature_squared + 1 + q->sigma * q->sigma)
    - 2 * q->etabar_squared_over_curvature_squared * (-q->spsi * q->torsion + q->I2 / q->B0) * q->G0 / q->B0;
}

/** Jacobian for the sigma equation.
 */
void Qsc::sigma_eq_jacobian(Vector& state, Matrix& jac, void* user_data) {
  // Get a pointer to the relevant Qsc object:
  Qsc* q = (Qsc*)user_data;
  
  // This function is always called after a call to the residual
  // function with the same state vector, so we don't need to extract
  // iota and sigma from the state vector here.
  int j;
  
  // d (Riccati equation) / d sigma:
  // For convenience we will fill all the columns now,
  // and re-write the first column in a moment.
  jac = q->d_d_varphi;
  for (j = 1; j < q->nphi; j++) {
    jac(j, j) += (q->iota + q->helicity * q->nfp) * 2 * q->sigma[j];    
  }

  // d (Riccati equation) / d iota:
  for (j = 0; j < q->nphi; j++) {
    jac(j, 0) = q->etabar_squared_over_curvature_squared[j]
      * q->etabar_squared_over_curvature_squared[j]
      + 1 + q->sigma[j] * q->sigma[j];
  }
}

void qsc::Qsc::solve_sigma_equation() {
  state = sigma0; // Initial guess for sigma
  state[0] = 0.0; // Initial guess for iota
    
  newton_solve(sigma_eq_residual, sigma_eq_jacobian,
	       state, residual, work1, work2, ipiv, work_matrix,
	       max_newton_iterations, max_linesearch_iterations,
	       newton_tolerance, verbose, this);

  if (verbose > 0) std::cout << "iota: " << iota << "  state[0]: " << state[0] << std::endl;
}
