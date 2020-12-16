#ifndef QSC_H
#define QSC_H

#ifndef QSC_REAL
#define QSC_REAL double
#endif

#include "vector_matrix.hpp"

namespace qsc {  

  const double pi = 3.141592653589793;
  
  Matrix differentiation_matrix(const int N, const QSC_REAL xmin, const QSC_REAL xmax);

  void hw();
  int return5();
}

#endif

