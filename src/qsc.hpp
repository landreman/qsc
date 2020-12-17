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
    Vector sinangle, cosangle, tempvec, tempvec1, tempvec2, tempvec3;
    Vector tangent_cylindrical1, tangent_cylindrical2, tangent_cylindrical3;
    Vector normal_cylindrical1, normal_cylindrical2, normal_cylindrical3;
    Vector binormal_cylindrical1, binormal_cylindrical2, binormal_cylindrical3;
    Vector d_tangent_d_l_cylindrical1, d_tangent_d_l_cylindrical2, d_tangent_d_l_cylindrical3;
    Vector torsion_numerator, torsion_denominator, B1Squared_over_curvatureSquared;
    
  public:
    Vector R0c, R0s, Z0c, Z0s;
    int nphi, nfp;
    double eta_bar;
    Vector phi, R0, Z0, R0p, Z0p, R0pp, Z0pp, R0ppp, Z0ppp;
    Vector d_l_d_phi, d2_l_d_phi2;
    Vector curvature, torsion;
    double d_phi, B0, G0, B0_over_abs_G0, abs_G0_over_B0, d_l_d_varphi;
    int sG, spsi;
    double axis_length, rms_curvature;
    double mean_R, mean_Z, standard_deviation_of_R, standard_deviation_of_Z;
    Matrix d_d_phi, d_d_varphi;
    Vector X1s, X1c;
    Vector Boozer_toroidal_angle;
    
    Qsc();
    void allocate();
    void init_axis();
  };
}

#endif

