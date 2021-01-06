#include <iostream>
#include <valarray>
#include "qsc.hpp"

using namespace qsc;

/** Use Newton's method to solve a nonlinear system of equations.
 *
 * The number of equations must equal the number of unknowns.
 *
 * The user is responsible for allocating all the Vectors and Matrix
 * used, so all these objects can be re-used for efficiency during a
 * parameter scan.
 *
 * @param state On entry, the initial guess for the solution. On exit,
 *        the solution.
 * @param residual Values on entry are ignored. On exit, the residuals
 *        at the solution.
 * @param state0 Merely a work array. Values on entry are ignored, and
 *        values on exit are irrelevant. Do not mistake this parameter
 *        for the initial guess.
 * @param step_direction Merely a work array. Values on entry are ignored, and
 *        values on exit are irrelevant.
 * @param ipiv Merely a work array. Values on entry are ignored, and
 *        values on exit are irrelevant.
 * @param m Merely a work array. Values on entry are ignored, and
 *        values on exit are irrelevant.
 * @param user_data This pointer allows you to pass any data you like to the
 *        residual and jacobian functions.
 */
int qsc::newton_solve(residual_function_type residual_function,
		       jacobian_function_type jacobian_function,
		       Vector& state,
		       Vector& residual,
		       Vector& state0, // Work array
		       Vector& step_direction, // Work array
		       std::valarray<int>& ipiv,
		       Matrix& m,
		       int max_newton_iterations,
		       int max_linesearch_iterations,
		       qscfloat tolerance,
		       int verbose,
		       void* user_data) {
  
  qscfloat tolerance_sq = tolerance * tolerance;

  // Find the initial residual:
  residual_function(state, residual, user_data);
  qscfloat initial_residual_norm_sq = dot_product(residual, residual);
  qscfloat residual_norm_sq = initial_residual_norm_sq;
  qscfloat last_residual_norm_sq, step_scale;
  if (verbose > 0) {
    std::cout << "Initial squared residual norm: " << initial_residual_norm_sq << std::endl;
  }

  int j_newton, j_linesearch;
  for (j_newton = 0; j_newton < max_newton_iterations; j_newton++) {
    last_residual_norm_sq = residual_norm_sq;
    if (residual_norm_sq < tolerance_sq) return NEWTON_CONVERGED;

    jacobian_function(state, m, user_data);

    state0 = state;
    if (verbose > 0) std::cout << "  Newton iteration " << j_newton << std::endl;
    // step_direction = - matrix \ residual
    // Note that LAPACK will
    // over-write step_direction with the solution, and over-write
    // the Jacobian with the LU factorization.
    step_direction = -residual;
    linear_solve(m, step_direction, ipiv);
    step_scale = 1.0;
    for (j_linesearch = 0; j_linesearch < max_linesearch_iterations; j_linesearch++) {
      state = state0 + step_scale * step_direction;
      residual_function(state, residual, user_data);
      residual_norm_sq = dot_product(residual, residual);
      if (verbose > 0) {
	std::cout << "    Line search step " << j_linesearch
		  << "  Squared residual norm: " << residual_norm_sq << std::endl;
      }
      if (residual_norm_sq < last_residual_norm_sq) break;
      step_scale /= 2.0;
    }
    if (residual_norm_sq > last_residual_norm_sq) {
      if (verbose > 0) std::cout << "Line search failed to reduce residual." << std::endl;
      // If the line search fails, stop the Newton iteration:
      state = state0;
      return NEWTON_LINESEARCH_FAILED;
      //break;
    }
  }
  
  if (residual_norm_sq < tolerance_sq) {
    return NEWTON_CONVERGED;
  } else {
    return NEWTON_MAX_ITERATIONS;
  }

}
