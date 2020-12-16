#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("2x2 differentiation matrix") {
  int n = 2;
  
  Matrix D1(n, n);
  D1 = 0.0;

  Matrix D2 = differentiation_matrix(n, 0, 2 * pi);

  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      CHECK(D1(k, j) == Approx(D2(k, j)));
    }
  }
}
