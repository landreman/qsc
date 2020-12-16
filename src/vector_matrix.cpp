#include <iostream>
#include <valarray>
#include <cassert>
#include "vector_matrix.hpp"

using namespace qsc;

Matrix::Matrix(index_type nrows_in, index_type ncols_in)
  : std::valarray<double>(nrows_in * ncols_in) // Call constructor of base class.
{
  nrows_ = nrows_in;
  ncols_ = ncols_in;
  len_ = nrows_ * ncols_;
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
  double sum = 0;
  for (int j = 0; j < m.nrows(); j++) {
    sum = 0;
    for (int k = 0; k < v.size(); k++) sum += m(j, k) * v[k];
    result[j] = sum;
  }
  return result;
}
