#include "doctest.h"
#include "qsc.hpp"
#include "vector_matrix.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("Create a matrix") {  
  Matrix m1(3, 4), m2(1, 3), m3(2, 2);
  CHECK(1 == 1);
}

TEST_CASE("Set matrix to a constant") {
  Matrix m1(3, 4);
  m1 = 42.0;
  for (int j = 0; j < m1.ncols(); j++) {
    for (int k = 0; k < m1.nrows(); k++) {
      CHECK(m1(k, j) == Approx(42.0));
    }
  }
}
/*
TEST_CASE("Set matrix to a constant") {
  Matrix m1(3, 4);
  m1 = 42.0;
  m1(1, 1) = 17.0;
  std::cout << m1(1, 1) << std::endl;
  
  std::cout << m1;

  m2 = 1.0;
  m3 = m2;
  // m3 = m1 + m2;

  Vector v1(4);
  v1 = -3.1;
  std::cout << "Here is v1: " << v1 << std::endl;

  Vector v2(4), v3(4);
  for (int j = 0; j < v2.size(); j++) v2[j] = j;
  //  v3 = (-v1) + 2 * v2;
  v3 = 2.0 * v2;
  std::cout << "Here is v3: " << v3 << std::endl;
}
*/

TEST_CASE("matrix-vector multiply 1") {
  Vector v1{0.6, 0.9, -1.2};
  Matrix m(2, 3);
  m(0, 0) = 1.1;
  m(1, 0) = -0.3;
  m(0, 1) = -0.1;
  m(1, 1) = 0.7;
  m(0, 2) = -1.3;
  m(1, 2) = 0.5;

  Vector v2 = m * v1;
  CHECK(v2[0] == Approx(2.13));
  CHECK(v2[1] == Approx(-0.15));
}

TEST_CASE("matrix-vector multiply 2") {
  Vector v1{0.6, 0.9};
  Matrix m(3, 2);
  m(0, 0) = 1.1;
  m(1, 0) = -0.3;
  m(2, 0) = -0.1;
  m(0, 1) = 0.7;
  m(1, 1) = -1.3;
  m(2, 1) = 0.5;

  Vector v2 = m * v1;
  CHECK(v2[0] == Approx(1.29));
  CHECK(v2[1] == Approx(-1.35));
  CHECK(v2[2] == Approx(0.39));
}

TEST_CASE("linear solve 1") {
  Vector v {0.6, -0.9};
  Matrix m(2, 2);
  m(0, 0) = 1.1;
  m(1, 0) = 0.7;
  m(0, 1) = -0.3;
  m(1, 1) = -1.3;

  std::valarray<int> IPIV(2);
  linear_solve(m, v, IPIV);
  
  CHECK(Approx(v[0]) == 0.860655737704918);
  CHECK(Approx(v[1]) == 1.15573770491803);
}

TEST_CASE("linear solve 3") {
  Vector v {0.6, -0.9, 0.4};
  Matrix m(3, 3);
  m(0, 0) = 1.1;
  m(1, 0) = 0.7;
  m(2, 0) = -0.2;
  m(0, 1) = -0.3;
  m(1, 1) = -1.3;
  m(2, 1) = 0.1;
  m(0, 2) = 0.6;
  m(1, 2) = -1.2;
  m(2, 2) = 2.1;

  std::valarray<int> IPIV(3);
  linear_solve(m, v, IPIV);
  
  CHECK(Approx(v[0]) == 0.661697247706422);
  CHECK(Approx(v[1]) == 0.852064220183486);
  CHECK(Approx(v[2]) == 0.212920489296636);
}
