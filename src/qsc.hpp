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

  class Qsc {
  private:
    Vector sinangle, cosangle;
  public:
    Vector R0c, R0s, Z0c, Z0s;
    int nphi, nfp;
    Vector phi, R0, Z0, R0p, Z0p, R0pp, Z0pp, R0ppp, Z0ppp;
    Vector curvature, torsion;
    double d_phi;

    void allocate();
    void init_axis();
  };
}

#endif

