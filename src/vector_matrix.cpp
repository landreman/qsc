#include <iostream>
#include <valarray>
#include <cassert>
#include "vector_matrix.hpp"

// Representation of BLAS and LAPACK routines we need:
extern "C" {
  // Matrix-vector multiply:
  void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA,
	      double* X, int* INCX, double* BETA, double* Y, int* INCY);

  void sgemv_(char* TRANS, int* M, int* N, float* ALPHA, float* A, int* LDA,
	      float* X, int* INCX, float* BETA, float* Y, int* INCY);

  // Matrix-matrix multiply:
  // void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA,
  //	      double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);

  // Solve linear system:
  // void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
}

// Choose either the single or double precision version of BLAS/LAPACK routines:
#ifdef SINGLE

#define gemv_ sgemv_

#else

#define gemv_ dgemv_

#endif

using namespace qsc;

Matrix::Matrix(index_type nrows_in, index_type ncols_in)
  : std::valarray<qscfloat>(nrows_in * ncols_in) // Call constructor of base class.
{
  nrows_ = nrows_in;
  ncols_ = ncols_in;
  len_ = nrows_ * ncols_;
}

void Matrix::resize(index_type nrows_in, index_type ncols_in, qscfloat v) {
  nrows_ = nrows_in;
  ncols_ = ncols_in;
  len_ = nrows_ * ncols_;
  std::valarray<qscfloat>::resize(nrows_ * ncols_, v);
}
	       
std::ostream& qsc::operator<< (std::ostream& os, Vector& v) {
  for (index_type j = 0; j < v.size(); j++) {
    if (j > 0) os << " ";
    os << v[j];
  }
  return os;
}

std::ostream& qsc::operator<< (std::ostream& os, Matrix& m) {
  for (index_type j = 0; j < m.nrows(); j++) {
    for (index_type k = 0; k < m.ncols(); k++) {
      if (k > 0) os << " ";
      os << m[j + m.nrows() * k];
    }
    os << std::endl;
  }
  return os;
}

Vector qsc::operator*(Matrix& m, Vector& v) {
  assert(m.ncols() == v.size());
  Vector result(m.nrows());

  /*
  // Non-BLAS simple-minded method:
  qscfloat sum = 0;
  for (int j = 0; j < m.nrows(); j++) {
    sum = 0;
    for (int k = 0; k < v.size(); k++) sum += m(j, k) * v[k];
    result[j] = sum;
  }
  */

  char TRANS = 'N';
  int INC = 1;
  int nrows = m.nrows();
  int ncols = m.ncols();
  qscfloat ALPHA = 1.0;
  qscfloat BETA = 0.0;
#ifdef SINGLE
  std::cout << "calling SINGLE precision BLAS" << std::endl;
#else
  std::cout << "calling DOUBLE precision BLAS" << std::endl;
#endif
  gemv_(&TRANS, &nrows, &ncols, &ALPHA, &m(0, 0), &nrows,
	      &v[0], &INC, &BETA, &result[0], &INC);
  return result;
}
