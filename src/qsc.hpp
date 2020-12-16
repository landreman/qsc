#ifndef QSC_H
#define QSC_H

#include <valarray>

#ifndef QSC_REAL
#define QSC_REAL double
#endif

namespace qsc {  

  const double pi = 3.141592653589793;
  
  typedef std::valarray<QSC_REAL> Vec;
  
  Vec differentiation_matrix(const int N, const QSC_REAL xmin, const QSC_REAL xmax);

  void hw();
  int return5();
}

#endif

