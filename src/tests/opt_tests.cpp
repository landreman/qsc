#include <stdexcept>
#include <iostream>
#include <gsl/gsl_vector.h>
#include "doctest.h"
#include "opt.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("Each term in the objective function should be approximately independent of nphi") {
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
  o1.weight_R0 = 4.0;
  o1.min_R0 = 0.8;
  o1.weight_d2_volume_d_psi2 = 5.0;
  o1.weight_XY2 = 6.0;
  o1.weight_XY2Prime = 7.0;
  o1.weight_XY3 = 8.0;
  o1.weight_XY3Prime = 9.0;
  o1.weight_grad_grad_B = 10.0;
  
  o2.weight_B20 = 2.0;
  o2.weight_iota = 3.0;
  o2.weight_R0 = 4.0;
  o2.min_R0 = 0.8;
  o2.weight_d2_volume_d_psi2 = 5.0;
  o2.weight_XY2 = 6.0;
  o2.weight_XY2Prime = 7.0;
  o2.weight_XY3 = 8.0;
  o2.weight_XY3Prime = 9.0;
  o2.weight_grad_grad_B = 10.0;

  o1.init();
  o2.init();
  
  gsl_vector *res1 = gsl_vector_alloc(o1.n_terms);
  gsl_vector *res2 = gsl_vector_alloc(o2.n_terms);
  o1.set_residuals(res1);
  o2.set_residuals(res2);

  qscfloat tol = 1.0e-8;
  CHECK(Approx(o1.objective_function).epsilon(tol) == o2.objective_function);
  CHECK(Approx(o1.B20_term).epsilon(tol) == o2.B20_term);
  CHECK(Approx(o1.iota_term).epsilon(tol) == o2.iota_term);
  CHECK(Approx(o1.R0_term).epsilon(1.0e-4) == o2.R0_term); // This term needs a wider tolerance
  CHECK(Approx(o1.d2_volume_d_psi2_term).epsilon(tol) == o2.d2_volume_d_psi2_term);
  CHECK(Approx(o1.XY2_term).epsilon(tol) == o2.XY2_term);
  CHECK(Approx(o1.XY2Prime_term).epsilon(tol) == o2.XY2Prime_term);
  CHECK(Approx(o1.XY3_term).epsilon(tol) == o2.XY3_term);
  CHECK(Approx(o1.XY3Prime_term).epsilon(tol) == o2.XY3Prime_term);
  CHECK(Approx(o1.grad_grad_B_term).epsilon(tol) == o2.grad_grad_B_term);

}
