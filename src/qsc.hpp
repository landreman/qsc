#ifndef QSC_H
#define QSC_H

#include <string>
#include "vector_matrix.hpp"

namespace qsc {  

  const qscfloat pi = 3.141592653589793;
  
  Matrix differentiation_matrix(const int N, const qscfloat xmin, const qscfloat xmax);

  typedef void (*residual_function_type)(Vector&, Vector&, void*);
  typedef void (*jacobian_function_type)(Vector&, Matrix&, void*);
  void newton_solve(residual_function_type, jacobian_function_type,
		    Vector&, Vector&, Vector&, Vector&, std::valarray<int>&,
		    Matrix&, int, int, qscfloat, int, void*);
  
  class Qsc {
  private:
    Vector sinangle, cosangle, tempvec, tempvec1, tempvec2, tempvec3;
    Vector tangent_cylindrical1, tangent_cylindrical2, tangent_cylindrical3;
    Vector normal_cylindrical1, normal_cylindrical2, normal_cylindrical3;
    Vector binormal_cylindrical1, binormal_cylindrical2, binormal_cylindrical3;
    Vector d_tangent_d_l_cylindrical1, d_tangent_d_l_cylindrical2, d_tangent_d_l_cylindrical3;
    Vector torsion_numerator, torsion_denominator, etabar_squared_over_curvature_squared;
    std::valarray<int> quadrant, ipiv;
    Vector state, residual, work1, work2;
    Matrix work_matrix;
    
    void calculate_helicity();
    static void sigma_eq_residual(Vector&, Vector&, void*);
    static void sigma_eq_jacobian(Vector&, Matrix&, void*);
    void calculate_grad_B_tensor();
    
  public:
    int verbose;
    Vector R0c, R0s, Z0c, Z0s;
    int nphi, nfp;
    qscfloat eta_bar, sigma0, B2c, B2s, I2, p2;
    Vector phi, R0, Z0, R0p, Z0p, R0pp, Z0pp, R0ppp, Z0ppp;
    Vector d_l_d_phi, d2_l_d_phi2;
    Vector curvature, torsion;
    qscfloat d_phi, B0, G0, B0_over_abs_G0, abs_G0_over_B0, d_l_d_varphi;
    int sG, spsi, helicity;
    qscfloat axis_length, rms_curvature;
    qscfloat mean_R, mean_Z, standard_deviation_of_R, standard_deviation_of_Z;
    Matrix d_d_phi, d_d_varphi;
    Vector X1s, X1c, sigma, Y1s, Y1c, elongation;
    Vector Boozer_toroidal_angle, L_grad_B;
    int max_newton_iterations, max_linesearch_iterations;
    qscfloat newton_tolerance;
    qscfloat iota, iota_N, grid_max_curvature, grid_max_elongation, mean_elongation;
    std::string order_r_option;
    Vector d_X1c_d_varphi, d_Y1s_d_varphi, d_Y1c_d_varphi;
    Rank3Tensor grad_B_tensor;
    qscfloat grid_min_L_grad_B;
    
    Qsc();
    Qsc(std::string);
    void input(std::string);
    void defaults();
    void allocate();
    void init_axis();
    void solve_sigma_equation();
    void r1_diagnostics();
    void calculate();
    void write_output();
  };
}

#endif

