#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;

TEST_CASE("2x2 differentiation matrix") {
  Vec D1 {
    0.0, 0.0,
    0.0, 0.0};

  Vec D2 = differentiation_matrix(2, 0, 2 * pi);
  
  // CHECK(D1 == D2);
}
