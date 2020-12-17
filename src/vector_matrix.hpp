#ifndef QSC_VECTOR_MATRIX_H
#define QSC_VECTOR_MATRIX_H

#include <valarray>
#include <iostream>

namespace qsc {
  
  typedef std::valarray<double> Vector;
  // typedef unsigned index_type;
  typedef int index_type;

  class Matrix : public std::valarray<double> {

  private:
    index_type nrows_, ncols_, len_;

  public:
    Matrix(index_type, index_type);
    index_type nrows();
    index_type ncols();
    // For info about matrix indexing:
    // https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
    double& operator()(index_type, index_type);
    double  operator()(index_type, index_type) const;
    Matrix& operator=(const double);
    Matrix& operator=(const Matrix&);
  };

  // These operators should be in the qsc namespace:
  // https://stackoverflow.com/questions/3891402/operator-overloading-and-namespaces
  // https://stackoverflow.com/questions/3623631/where-should-non-member-operator-overloads-be-placed
  std::ostream& operator<<(std::ostream&, Vector&);
  std::ostream& operator<<(std::ostream&, Matrix&);
  Vector operator*(Matrix&, Vector&);

  /*
  // Multiplication of an int with std::valarray<double> is not included in some compilers:
  inline Vector operator*(int j, Vector& v) {
    return std::operator*(double(j), v);
  }
  */

  // inline functions must be included in every file that uses them,
  // so these functions should go in the header file.
  inline index_type Matrix::nrows() {
    return nrows_;
  }

  inline index_type Matrix::ncols() {
    return ncols_;
  }

  inline double& Matrix::operator()(index_type m, index_type n) {
    return (*this)[m + nrows_ * n];
  }

  inline double Matrix::operator()(index_type m, index_type n) const {
    return (*this)[m + nrows_ * n];
  }

  inline Matrix& Matrix::operator=(const double v) {
    // Delegate to parent class:
    std::valarray<double>::operator=(v);
    return *this;
  }

  inline Matrix& Matrix::operator=(const Matrix &v) {
    // Delegate to parent class:
    std::valarray<double>::operator=(v);
    return *this;
  }


} // namespace qsc


#endif
