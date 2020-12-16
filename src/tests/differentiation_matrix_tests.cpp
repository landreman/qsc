#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("differentiation matrix: diag and row/col sums") {
  double sum;
  for (int n = 2; n < 10; n++) {
    for (int j_min = 0; j_min < 5; j_min++) {
      double xmin = (j_min - 3.5) * 0.7; 
      for (int j_max = 0; j_max < 4; j_max++) {
	double xmax = xmin + (j_max + 1) * 0.6;
	Matrix ddx = differentiation_matrix(n, xmin, xmax);

	// Diagonal should be 0:
	for (int j = 0; j < n; j++) {
	  CHECK(Approx(ddx(j, j)) == 0.0);
	}

	// Row sums should be 0:
	for (int j = 0; j < n; j++) {
	  sum = 0;
	  for (int k = 0; k < n; k++) sum += ddx(j, k);
	  CHECK(Approx(sum) == 0.0);
	}

	// Column sums should be 0:
	for (int j = 0; j < n; j++) {
	  sum = 0;
	  for (int k = 0; k < n; k++) sum += ddx(k, j);
	  CHECK(Approx(sum) == 0.0);
	}
	
      }
    }
  }
}

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

TEST_CASE("3x3 differentiation matrix") {
  int n = 3;
  
  Matrix D1(n, n);
  D1 = 0.0;
  double x = 0.577350269189626;
  D1(1, 0) = -x;
  D1(2, 1) = -x;
  D1(0, 2) = -x;
  D1(0, 1) = x;
  D1(1, 2) = x;
  D1(2, 0) = x;

  Matrix D2 = differentiation_matrix(n, 0, 2 * pi);

  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      CHECK(D1(k, j) == Approx(D2(k, j)));
    }
  }
}

TEST_CASE("4x4 differentiation matrix") {
  int n = 4;
  
  Matrix D1(n, n);
  D1 = 0.0;
  D1(1, 0) = -0.5;
  D1(2, 1) = -0.5;
  D1(3, 2) = -0.5;
  D1(0, 3) = -0.5;
  D1(0, 1) = 0.5;
  D1(1, 2) = 0.5;
  D1(2, 3) = 0.5;
  D1(3, 0) = 0.5;

  Matrix D2 = differentiation_matrix(n, 0, 2 * pi);
  
  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      CHECK(D1(k, j) == Approx(D2(k, j)));
    }
  }
}

TEST_CASE("5x5 differentiation matrix") {
  int n = 5;
  
  Matrix D1(n, n);
  D1 = 0.0;
  double e = 0.85065080835204;
  double f = 0.525731112119134;
  for (int j = 0; j < n; j++) {
    D1(j, (j + 1) % n) = e;
    D1(j, (j + 2) % n) = -f;
    D1(j, (j + 3) % n) = f;
    D1(j, (j + 4) % n) = -e;
  }
  
  Matrix D2 = differentiation_matrix(n, 0, 2 * pi);

  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      CHECK(D1(k, j) == Approx(D2(k, j)));
    }
  }
}

TEST_CASE("2x2 differentiation matrix, shifted") {
  int n = 2;
  
  Matrix D1(n, n);
  D1 = 0.0;

  Matrix D2 = differentiation_matrix(n, -2.1, 3.7);

  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      CHECK(D1(k, j) == Approx(D2(k, j)));
    }
  }
}

TEST_CASE("3x3 differentiation matrix, shifted") {
  int n = 3;
  
  Matrix D1(n, n);
  D1 = 0.0;
  double x = 0.625448056632489;
  D1(1, 0) = -x;
  D1(2, 1) = -x;
  D1(0, 2) = -x;
  D1(0, 1) = x;
  D1(1, 2) = x;
  D1(2, 0) = x;

  Matrix D2 = differentiation_matrix(n, -2.1, 3.7);

  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      CHECK(D1(k, j) == Approx(D2(k, j)));
    }
  }
}

TEST_CASE("4x4 differentiation matrix, shifted") {
  int n = 4;
  
  Matrix D1(n, n);
  D1 = 0.0;
  double x = 0.541653905791344;
  D1(1, 0) = -x;
  D1(2, 1) = -x;
  D1(3, 2) = -x;
  D1(0, 3) = -x;
  D1(0, 1) = x;
  D1(1, 2) = x;
  D1(2, 3) = x;
  D1(3, 0) = x;

  Matrix D2 = differentiation_matrix(n, -2.1, 3.7);
  
  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      CHECK(D1(k, j) == Approx(D2(k, j)));
    }
  }
}

TEST_CASE("5x5 differentiation matrix, shifted") {
  int n = 5;
  
  Matrix D1(n, n);
  D1 = 0.0;
  double e = 0.921516665616892;
  double f = 0.569528620550711;
  for (int j = 0; j < n; j++) {
    D1(j, (j + 1) % n) = e;
    D1(j, (j + 2) % n) = -f;
    D1(j, (j + 3) % n) = f;
    D1(j, (j + 4) % n) = -e;
  }
  
  Matrix D2 = differentiation_matrix(n, -2.1, 3.7);

  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      CHECK(D1(k, j) == Approx(D2(k, j)));
    }
  }
}
