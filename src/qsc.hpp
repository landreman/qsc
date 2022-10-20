#ifndef QSC_H
#define QSC_H

#include <string>
#include "vector_matrix.hpp"

namespace qsc {  

  const qscfloat pi = 3.141592653589793;
  const qscfloat mu0 = (4.0e-7) * pi;
  
  typedef unsigned long long int big;
  
  Matrix differentiation_matrix(const int N, const qscfloat xmin, const qscfloat xmax);

  typedef void (*residual_function_type)(Vector&, Vector&, void*);
  typedef void (*jacobian_function_type)(Vector&, Matrix&, void*);
  int newton_solve(residual_function_type, jacobian_function_type,
		    Vector&, Vector&, Vector&, Vector&, std::valarray<int>&,
		    Matrix&, int, int, qscfloat, int, void*);

  const std::string ORDER_R_OPTION_R1 = "r1";
  const std::string ORDER_R_OPTION_R2 = "r2";
  const std::string ORDER_R_OPTION_R2p1 = "r2.1";

  int driver(int, char**);

  enum {
    NEWTON_CONVERGED,
    NEWTON_MAX_ITERATIONS,
    NEWTON_LINESEARCH_FAILED};
  
  class Qsc {
  private:
    Vector sinangle, cosangle, tempvec, tempvec1, tempvec2, tempvec3;
    Vector tangent_cylindrical1, tangent_cylindrical2, tangent_cylindrical3;
    Vector normal_cylindrical1, normal_cylindrical2, normal_cylindrical3;
    Vector binormal_cylindrical1, binormal_cylindrical2, binormal_cylindrical3;
    Vector d_tangent_d_l_cylindrical1, d_tangent_d_l_cylindrical2, d_tangent_d_l_cylindrical3;
    Vector torsion_numerator, torsion_denominator, etabar_squared_over_curvature_squared;
    std::valarray<int> quadrant, ipiv, r2_ipiv;
    Vector state, residual, work1, work2;
    Matrix work_matrix;
    Vector V1, V2, V3, rc, rs, qc, qs, r2_rhs;
    Matrix r2_matrix;
    Vector Y2s_from_X20, Y2s_inhomogeneous, Y2c_from_X20, Y2c_inhomogeneous;
    Vector fX0_from_X20, fX0_from_Y20, fX0_inhomogeneous;
    Vector fXs_from_X20, fXs_from_Y20, fXs_inhomogeneous;
    Vector fXc_from_X20, fXc_from_Y20, fXc_inhomogeneous;
    Vector fY0_from_X20, fY0_from_Y20, fY0_inhomogeneous;
    Vector fYs_from_X20, fYs_from_Y20, fYs_inhomogeneous;
    Vector fYc_from_X20, fYc_from_Y20, fYc_inhomogeneous;
    
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
    Vector Boozer_toroidal_angle, L_grad_B, L_grad_B_inverse;
    Vector L_grad_grad_B, L_grad_grad_B_inverse;
    int max_newton_iterations, max_linesearch_iterations;
    qscfloat newton_tolerance, grid_min_R0, G2, I2_over_B0;
    qscfloat iota, iota_N, grid_max_curvature, grid_max_elongation, mean_elongation;
    std::string order_r_option;
    bool at_least_order_r2, order_r2p1, order_r3;
    Vector X20, X2s, X2c, Y20, Y2s, Y2c, Z20, Z2s, Z2c;
    Vector d_X1c_d_varphi, d_Y1s_d_varphi, d_Y1c_d_varphi, B20, B20_anomaly;
    Rank3Tensor grad_B_tensor;
    Rank4Tensor grad_grad_B_tensor;
    qscfloat grid_min_L_grad_B, grid_min_L_grad_grad_B, beta_1s, B20_mean, B20_residual, B20_grid_variation;
    Vector d_curvature_d_varphi, d_torsion_d_varphi;
    Vector d_X20_d_varphi, d_X2s_d_varphi, d_X2c_d_varphi;
    Vector d_Y20_d_varphi, d_Y2s_d_varphi, d_Y2c_d_varphi;
    Vector d_Z20_d_varphi, d_Z2s_d_varphi, d_Z2c_d_varphi;
    Vector d2_X1c_d_varphi2, d2_Y1c_d_varphi2, d2_Y1s_d_varphi2;
    Vector d2_X20_d_varphi2, d2_X2s_d_varphi2, d2_X2c_d_varphi2;
    Vector d2_Y20_d_varphi2, d2_Y2s_d_varphi2, d2_Y2c_d_varphi2;
    Vector d2_Z20_d_varphi2, d2_Z2s_d_varphi2, d2_Z2c_d_varphi2;
    qscfloat d2_volume_d_psi2, DGeod_times_r2, DWell_times_r2, DMerc_times_r2;
    qscfloat r_singularity_robust, beta;
    Vector r_hat_singularity_robust;
    int newton_result;
    Vector X3c1, X3c3, X3s1, X3s3, Y3c1, Y3c3, Y3s1, Y3s3;
    Vector Z3c1, Z3c3, Z3s1, Z3s3, lambda_for_XY3;
    Vector d_X3c1_d_varphi, d_X3c3_d_varphi, d_X3s1_d_varphi, d_X3s3_d_varphi;
    Vector d_Y3c1_d_varphi, d_Y3c3_d_varphi, d_Y3s1_d_varphi, d_Y3s3_d_varphi;
    Vector d2_X3c1_d_varphi2, d2_X3c3_d_varphi2, d2_X3s1_d_varphi2, d2_X3s3_d_varphi2;
    Vector d2_Y3c1_d_varphi2, d2_Y3c3_d_varphi2, d2_Y3s1_d_varphi2, d2_Y3s3_d_varphi2;
    qscfloat grid_max_XY2, grid_max_Z2, grid_max_XY3;
    qscfloat grid_max_d_XY2_d_varphi, grid_max_d_Z2_d_varphi, grid_max_d_XY3_d_varphi;
    qscfloat grid_max_d2_XY2_d_varphi2, grid_max_d2_XY3_d_varphi2;
    
    Qsc();
    Qsc(std::string);
    void input(std::string);
    void defaults();
    void validate();
    void allocate();
    void init_axis();
    void solve_sigma_equation();
    void r1_diagnostics();
    void calculate_r2();
    void calculate_r2p1();
    void r2_diagnostics();
    void init();
    void calculate();
    Rank4Tensor calculate_grad_grad_B_tensor_alt();
    void write_netcdf(std::string);
    void read_netcdf(std::string, char);
    void run(std::string);
    void calculate_grad_grad_B_tensor();
    void mercier();
    void calculate_r_singularity();
  };
  
  std::string outfile(std::string);
}

#endif
