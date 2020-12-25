#include <iostream>
#include <valarray>
#include <ctime>
#include <chrono>
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

void Qsc::solve_sigma_equation() {
  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }
  
  state = sigma0; // Initial guess for sigma
  state[0] = 0.0; // Initial guess for iota
    
  newton_solve(sigma_eq_residual, sigma_eq_jacobian,
	       state, residual, work1, work2, ipiv, work_matrix,
	       max_newton_iterations, max_linesearch_iterations,
	       newton_tolerance, verbose, this);

  if (verbose > 0) {
    end_time = std::clock();
    auto end = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for sigma equation from chrono:           "
              << elapsed.count() << " seconds" << std::endl;
    std::cout << "Time for sigma equation from ctime (CPU time): "
              << double(end_time - start_time) / CLOCKS_PER_SEC
              << " seconds" << std::endl;
  }

}

void Qsc::r1_diagnostics() {
  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }

  iota_N = iota + helicity * nfp;
  I2_over_B0 = I2 / B0;
  
  Y1s = (sG * spsi / eta_bar) * curvature;
  Y1c = Y1s * sigma;

  // Use (R,Z) for elongation in the (R,Z) plane
  // or use (X,Y) for elongation in the plane perpendicular to the magnetic axis.
  // tempvec1 = p, tempvec2 = q
  tempvec1 = X1s * X1s + X1c * X1c + Y1s * Y1s + Y1c * Y1c;
  tempvec2 = X1s * Y1c - X1c * Y1s;
  elongation = (tempvec1 + sqrt(tempvec1 * tempvec1 - 4 * tempvec2 * tempvec2))
    / (2 * abs(tempvec2));

  grid_max_elongation = elongation.max();
  tempvec = elongation * d_l_d_phi;
  mean_elongation = tempvec.sum() / d_l_d_phi.sum();
  //index = np.argmax(elongation);
  //max_elongation = -fourier_minimum(-elongation);

  matrix_vector_product(d_d_varphi, X1c, d_X1c_d_varphi);
  matrix_vector_product(d_d_varphi, Y1s, d_Y1s_d_varphi);
  matrix_vector_product(d_d_varphi, Y1c, d_Y1c_d_varphi);

  calculate_grad_B_tensor();

  if (verbose > 0) {
    end_time = std::clock();
    auto end = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for r1_diagnostics from chrono:           "
              << elapsed.count() << " seconds" << std::endl;
    std::cout << "Time for r1_diagnostics from ctime (CPU time): "
              << double(end_time - start_time) / CLOCKS_PER_SEC
              << " seconds" << std::endl;
  }
}
