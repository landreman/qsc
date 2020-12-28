#include "doctest.h"
#include "quartic_roots.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("Quartic roots 1") {
  qscfloat coefficients[] = {0.2, -0.1, -0.8, 0.5, -0.7};
  qscfloat real_parts[4], imag_parts[4];
  
  quartic_roots(coefficients, real_parts, imag_parts);
  /*
  std::cout << "Roots found:" << std::endl;
  for (int j = 0; j < 4; j++) {
    std::cout << "  " << real_parts[j] << " + " << imag_parts[j] << "i" << std::endl;
  }
  */

  // I'm not sure the order of the roots returned is the same for all
  // implementations of LAPACK. This complicates checking the answer.
  int root1found = 0, root2found = 0, root3found = 0, root4found = 0;
  qscfloat tol;
  if (single) {
    tol = 1.0e-6;
  } else {
    tol = 1.0e-12;
  }
  for (int j = 0; j < 4; j++) {
    if(std::abs(-2.1867382836778 - real_parts[j]) < tol &&
       std::abs(0.0              - imag_parts[j]) < tol) root1found++;

    if(std::abs(2.1617999494632 - real_parts[j]) < tol &&
       std::abs(0.0             - imag_parts[j]) < tol) root2found++;

    if(std::abs(0.262469167107297 - real_parts[j]) < tol &&
       std::abs(0.819445992779646 - imag_parts[j]) < tol) root3found++;

    if(std::abs( 0.262469167107297 - real_parts[j]) < tol &&
       std::abs(-0.819445992779646 - imag_parts[j]) < tol) root4found++;
  }
  CHECK(root1found == 1);
  CHECK(root2found == 1);
  CHECK(root3found == 1);
  CHECK(root4found == 1);
}

TEST_CASE("Quartic roots 2") {
  qscfloat coefficients[] = {1.0, 0.0, -7.0, -6.0, 0.0};
  qscfloat real_parts[4], imag_parts[4];
  
  quartic_roots(coefficients, real_parts, imag_parts);
  /*
  std::cout << "Roots found:" << std::endl;
  for (int j = 0; j < 4; j++) {
    std::cout << "  " << real_parts[j] << " + " << imag_parts[j] << "i" << std::endl;
  }
  */

  // I'm not sure the order of the roots returned is the same for all
  // implementations of LAPACK. This complicates checking the answer.
  int root1found = 0, root2found = 0, root3found = 0, root4found = 0;
  qscfloat tol;
  if (single) {
    tol = 1.0e-6;
  } else {
    tol = 1.0e-12;
  }
  for (int j = 0; j < 4; j++) {
    if(std::abs(0.0 - real_parts[j]) < tol &&
       std::abs(0.0 - imag_parts[j]) < tol) root1found++;

    if(std::abs(-1.0 - real_parts[j]) < tol &&
       std::abs( 0.0 - imag_parts[j]) < tol) root2found++;

    if(std::abs(-2.0 - real_parts[j]) < tol &&
       std::abs( 0.0 - imag_parts[j]) < tol) root3found++;

    if(std::abs( 3.0 - real_parts[j]) < tol &&
       std::abs( 0.0 - imag_parts[j]) < tol) root4found++;
  }
  CHECK(root1found == 1);
  CHECK(root2found == 1);
  CHECK(root3found == 1);
  CHECK(root4found == 1);
}
