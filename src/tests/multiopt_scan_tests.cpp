#include <iostream>
#include "doctest.h"
#include "multiopt_scan.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("Check linear and logarithmic spacing of scan parameters. [multiopt_scan]") {
  if (single) return;
  MultiOptScan mos;
  mos.params = {"R0c1", "weight_B20"};
  mos.params_min = {0.1, 0.01};
  mos.params_max = {0.4, 1.0};
  mos.params_n = {4, 3};
  mos.params_log = {false, true};

  mos.init();
  CHECK(mos.n_scan_all == 12);
  
  CHECK(Approx(mos.params_vals[0][0]) == 0.1);
  CHECK(Approx(mos.params_vals[0][1]) == 0.2);
  CHECK(Approx(mos.params_vals[0][2]) == 0.3);
  CHECK(Approx(mos.params_vals[0][3]) == 0.4);
  
  CHECK(Approx(mos.params_vals[1][0]) == 0.01);
  CHECK(Approx(mos.params_vals[1][1]) == 0.1);
  CHECK(Approx(mos.params_vals[1][2]) == 1.0);
}
