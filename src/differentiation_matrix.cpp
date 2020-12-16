#include <cmath>
#include "qsc.hpp"

using namespace qsc;

/**
 * Return the spectral differentiation matrix for n grid points on the
 * periodic domain [xmax, xmax). This routine is based on the matlab
 * code in the DMSuite package by S.C. Reddy and J.A.C. Weideman,
 * available at
 * http://www.mathworks.com/matlabcentral/fileexchange/29
 * or here:
 * http://dip.sun.ac.za/~weideman/research/differ.html  
 */
Matrix qsc::differentiation_matrix(const int N, const QSC_REAL xmin, const QSC_REAL xmax) {
    
  double h = 2 * pi / N;

  int n1 = floor((N - 1.0) / 2);
  int n2 = ceil((N - 1.0) / 2);
  Vector topc(n2);
  Vector col1(N);
  Matrix ddx(N, N);
  col1[0] = 0.0;
  int j, k;
    
  if (N % 2 == 0) {
    // Even size:
    
    for (j = 0; j < n2; j++) {
      topc[j] = 0.5 / tan((j + 1) * h / 2);
    }
    for (j = 1; j < n2 + 1; j++) {
      col1[j] = topc[j - 1];
    }
    for (j = 0; j < n1; j++) {
      col1[n2 + 1 + j] = -topc[n1 - 1 - j];
    }
    
  } else {
    // Odd size:
    
    for (j = 0; j < n2; j++) {
      topc[j] = 0.5 / sin((j + 1) * h / 2);
    }
    for (j = 1; j < n2 + 1; j++) {
      col1[j] = topc[j - 1];
    }
    for (j = 0; j < n1; j++) {
      col1[n2 + 1 + j] = topc[n1 - 1 - j];
    }
  }
  
  for (j = 1; j < N; j += 2) {
    col1[j] = -col1[j];
  }
  col1 = 2 * pi / (xmax - xmin) * col1;

  // Create a Toeplitz matrix:
  for (j = 0; j < N; j++) {
    for (k = j; k < N; k++) {
      ddx(j, k) = -col1[k - j];
    }
    for (k = 0; k < j; k++) {
      ddx(j, k) = col1[j - k];
    }
  }
  
  return ddx;
}
