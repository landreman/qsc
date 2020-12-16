#include "doctest.h"
#include "qsc.hpp"
#include "vector_matrix.hpp"

using namespace qsc;

TEST_CASE("Create a matrix") {  
  Matrix m1(3, 4), m2(1, 3), m3(2, 2);
  CHECK(1 == 1);
}

TEST_CASE("Set matrix to a constant") {
  Matrix m1(3, 4);
  m1 = 42.0;
  for (int j = 0; j < m1.ncols(); j++) {
    for (int k = 0; k < m1.nrows(); k++) {
      CHECK(m1(k, j) == 42.0);
    }
  }
  //CHECK(1 == 1);
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
