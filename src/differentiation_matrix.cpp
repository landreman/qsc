#include <cmath>
#include "qsc.hpp"

namespace qsc {
  
  Vec differentiation_matrix(const int N, const QSC_REAL xmin, const QSC_REAL xmax) {
    
    double h = 2 * pi / N;

    // Do floor and ceil need a "std::" in front?
    int n1 = floor((N - 1.0) / 2);
    int n2 = ceil((N - 1.0) / 2);
    Vec topc(n2);
    Vec col1(N);
    Vec ddx(N * N);
    col1[0] = 0.0;
    int j, k;
    
    if (N % 2 == 0) {
      for (j = 0; j < n2; j++) {
	topc[j] = 0.5 / tan((j + 1) * h / 2);
      }
      for (j = 1; j < n2 + 1; j++) {
	col1[j] = topc[j - 1];
      }
      for (j = 0; j < n1; j++) {
	col1[n2 + 1 + j] = topc[n1 - 1 - j];
      }
      for (j = 1; j < N; j += 2) {
	col1[j] = -col1[j];
      }
      col1 = 2 * pi / (xmax - xmin) * col1;
      // Create a Toeplitz matrix:
      for (j = 0; j < N; j++) {
	for (k = j; k < N; k++) {
	  ddx[j + N * k] = -col1[k - j];
	}
	for (k = 0; k < j; k++) {
	  ddx[j + N * k] = col1[j - k];
	}
      }
    } else {
      for (j = 0; j < n2; j++) {
	topc[j] = 0.5 / sin((j + 1) * h / 2);
      }
    }
    
    return ddx;
  }
  
}
