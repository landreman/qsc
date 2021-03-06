#ifndef QSC_VECTOR_MATRIX_H
#define QSC_VECTOR_MATRIX_H

#include <valarray>
#include <iostream>

namespace qsc {

#ifdef SINGLE
  typedef float qscfloat;
  const int single = 1;
#else
  typedef double qscfloat;
  const int single = 0;
#endif
  
  typedef std::valarray<qscfloat> Vector;
  // typedef unsigned index_type;
  typedef int index_type;

  class Matrix : public std::valarray<qscfloat> {

  private:
    index_type nrows_, ncols_, len_;

  public:
    Matrix();
    Matrix(index_type, index_type);
    index_type nrows();
    index_type ncols();
    void resize(index_type, index_type, qscfloat);
    // For info about matrix indexing:
    // https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
    qscfloat& operator()(index_type, index_type);
    qscfloat  operator()(index_type, index_type) const;
    Matrix& operator=(const qscfloat);
    Matrix& operator=(const Matrix&);
  };

  void matrix_vector_product(Matrix&, Vector&, Vector&);
  qscfloat dot_product(Vector&, Vector&);
  void linear_solve(Matrix&, Vector&, std::valarray<int>&);
  
  // These operators should be in the qsc namespace:
  // https://stackoverflow.com/questions/3891402/operator-overloading-and-namespaces
  // https://stackoverflow.com/questions/3623631/where-should-non-member-operator-overloads-be-placed
  std::ostream& operator<<(std::ostream&, Vector&);
  std::ostream& operator<<(std::ostream&, Matrix&);
  std::ostream& operator<<(std::ostream&, std::valarray<bool>&);
  Vector operator*(Matrix&, Vector&);

  
  // Multiplication of an int with std::valarray<qscfloat> is not included in some compilers:
  inline Vector operator*(int j, Vector& v) {
    //return qscfloat(j) * v;
    return std::operator*(qscfloat(j), v);
  }

  // gcc complains about ambiguity if we do not have this next little function:
  inline Vector operator*(qscfloat j, Vector& v) {
    return std::operator*(j, v);
  }
  
  /*
  inline Vector operator*(float j, Vector& v) {
    return std::operator*(qscfloat(j), v);
  }
  
  inline Vector operator*(double j, Vector& v) {
    return std::operator*(qscfloat(j), v);
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

  inline qscfloat& Matrix::operator()(index_type m, index_type n) {
    return (*this)[m + nrows_ * n];
  }

  inline qscfloat Matrix::operator()(index_type m, index_type n) const {
    return (*this)[m + nrows_ * n];
  }

  inline Matrix& Matrix::operator=(const qscfloat v) {
    // Delegate to parent class:
    std::valarray<qscfloat>::operator=(v);
    return *this;
  }

  inline Matrix& Matrix::operator=(const Matrix &v) {
    // Delegate to parent class:
    std::valarray<qscfloat>::operator=(v);
    return *this;
  }

  //////////////////////////////////////////////////////////
  
  class Rank3Tensor : public std::valarray<qscfloat> {

  private:
    index_type d1_, d2_, d3_, len_;

  public:
    Rank3Tensor();
    Rank3Tensor(index_type, index_type, index_type);
    void resize(index_type, index_type, index_type, qscfloat);
    // For info about matrix indexing:
    // https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
    qscfloat& operator()(index_type, index_type, index_type);
    qscfloat  operator()(index_type, index_type, index_type) const;
    void set_row(Vector&, index_type, index_type);
  };

  inline qscfloat& Rank3Tensor::operator()(index_type j1, index_type j2, index_type j3) {
    return (*this)[j1 + d1_ * (j2 + d2_ * j3)];
  }

  inline qscfloat Rank3Tensor::operator()(index_type j1, index_type j2, index_type j3) const {
    return (*this)[j1 + d1_ * (j2 + d2_ * j3)];
  }

  //////////////////////////////////////////////////////////
  
  class Rank4Tensor : public std::valarray<qscfloat> {

  private:
    index_type d1_, d2_, d3_, d4_, len_;

  public:
    Rank4Tensor();
    Rank4Tensor(index_type, index_type, index_type, index_type);
    void resize(index_type, index_type, index_type, index_type, qscfloat);
    // For info about matrix indexing:
    // https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
    qscfloat& operator()(index_type, index_type, index_type, index_type);
    qscfloat  operator()(index_type, index_type, index_type, index_type) const;
    void set_row(Vector&, index_type, index_type, index_type);
  };

  inline qscfloat& Rank4Tensor::operator()(index_type j1, index_type j2, index_type j3, index_type j4) {
    return (*this)[j1 + d1_ * (j2 + d2_ * (j3 + d3_ * j4))];
  }

  inline qscfloat Rank4Tensor::operator()(index_type j1, index_type j2, index_type j3, index_type j4) const {
    return (*this)[j1 + d1_ * (j2 + d2_ * (j3 + d3_ * j4))];
  }


} // namespace qsc


#endif
