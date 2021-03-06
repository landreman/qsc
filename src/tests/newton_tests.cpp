#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

void residual_function1(Vector& state, Vector& residual, void* user_data) {
  residual[0] = state[1] - exp(state[0]);
  residual[1] = state[0] + state[1];
}

void jacobian_function1(Vector& state, Matrix& jac, void* user_data) {
  jac(0, 0) = -exp(state[0]);
  jac(0, 1) = 1.0;
  jac(1, 0) = 1.0;
  jac(1, 1) = 1.0;
}

TEST_CASE("Newton solve") {
  Matrix m(2, 2);
  Vector state {0.0, 0.0};
  Vector work1(2), work2(2), residual(2);
  std::valarray<int> ipiv(2);

  int max_newton_iterations = 10;
  int max_linesearch_iterations = 5;
  qscfloat tolerance = 1.0e-12;
  if (single) tolerance = 1.0e-7;
  int verbose = 1;
  int result;

  result = newton_solve(residual_function1, jacobian_function1,
			state, residual, work1, work2, ipiv, m,
			max_newton_iterations, max_linesearch_iterations,
			tolerance, verbose, NULL);
  CHECK(Approx(state[0]) == -0.5671432904097838);
  CHECK(Approx(state[1]) == 0.5671432904097838);
  CHECK(result == NEWTON_CONVERGED);
}

