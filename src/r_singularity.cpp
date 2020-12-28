#include <ctime>
#include <chrono>
#include <cmath>
#include "qsc.hpp"
#include "quartic_roots.hpp"

using namespace qsc;

/** Compute \hat{r}_c(\varphi) from section 4 of Landreman, J Plasma Physics (2021).
 */
void Qsc::calculate_r_singularity() {
  qscfloat lp = abs_G0_over_B0; // shorthand
  int j;
  qscfloat K0, K2s, K2c, K4s, K4c;
  qscfloat coefficients[5], real_parts[4], imag_parts[4];
  qscfloat g0, g1c, g20, g2s, g2c;
  qscfloat g3s1, g3s3, g3c1, g3c3, g40, g4s2, g4s4, g4c2, g4c4;
  qscfloat rc, sin2theta, abs_cos2theta, residual_if_varpi_plus, residual_if_varpi_minus, cos2theta;
  int varpi, jr;
  qscfloat abs_costheta, abs_sintheta, costheta, sintheta, sintheta_at_rc, costheta_at_rc;
  qscfloat quadratic_A, quadratic_B, quadratic_C, radical, rr, residual;
  int varsigma, sign_quadratic;
  qscfloat sin2_cos2_1_tol, acceptable_residual;
  bool get_cos_from_cos2;

  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }

  if (single) {
    sin2_cos2_1_tol = 1.0e-6;
    acceptable_residual = 3.0e-2; // 1e-3
  } else {
    // Double precision
    sin2_cos2_1_tol = 1.0e-13;
    acceptable_residual = 1.0e-5;
  }
  
  for (j = 0; j < nphi; j++) {
    if (verbose > 1) std::cout << "---- r_singularity calculation for jphi = " << j << " ----" << std::endl;
    
    // Write sqrt(g) = r * [g0 + r*g1c*cos(theta) + (r^2)*(g20 + g2s*sin(2*theta) + g2c*cos(2*theta) + ...]
    // The coefficients are evaluated in "20200322-02 Max r for Garren Boozer.nb", in the section "Order r^2 construction, quasisymmetry"

    g0 = lp * X1c[j] * Y1s[j];

    // g1s = -2*X20[j]*Y1c[j] + 2*X2c[j]*Y1c[j] + 2*X2s[j]*Y1s[j] + 2*X1c[j]*Y20[j] - 2*X1c[j]*Y2c[j]
    // g1s vanishes for quasisymmetry.

    g1c = lp*(-2*X2s[j]*Y1c[j] + 2*X20[j]*Y1s[j] + 2*X2c[j]*Y1s[j] + 2*X1c[j]*Y2s[j] - X1c[j]*X1c[j]*Y1s[j]*curvature[j]);

    g20 = -4*lp*X2s[j]*Y2c[j] + 4*lp*X2c[j]*Y2s[j] + lp*X1c[j]*X2s[j]*Y1c[j]*curvature[j] - 
      2*lp*X1c[j]*X20[j]*Y1s[j]*curvature[j] - lp*X1c[j]*X2c[j]*Y1s[j]*curvature[j] - 
      lp*X1c[j]*X1c[j]*Y2s[j]*curvature[j] + 2*lp*Y1c[j]*Y1s[j]*Z2c[j]*torsion[j] - 
      lp*X1c[j]*X1c[j]*Z2s[j]*torsion[j] - lp*Y1c[j]*Y1c[j]*Z2s[j]*torsion[j] + lp*Y1s[j]*Y1s[j]*Z2s[j]*torsion[j] - 
      Y1s[j]*Z20[j]*d_X1c_d_varphi[j] - Y1s[j]*Z2c[j]*d_X1c_d_varphi[j] + 
      Y1c[j]*Z2s[j]*d_X1c_d_varphi[j] - X1c[j]*Z2s[j]*d_Y1c_d_varphi[j] - 
      X1c[j]*Z20[j]*d_Y1s_d_varphi[j] + X1c[j]*Z2c[j]*d_Y1s_d_varphi[j] + 
      X1c[j]*Y1s[j]*d_Z20_d_varphi[j];

    g2c = -4*lp*X2s[j]*Y20[j] + 4*lp*X20[j]*Y2s[j] + 
      lp*X1c[j]*X2s[j]*Y1c[j]*curvature[j] - lp*X1c[j]*X20[j]*Y1s[j]*curvature[j] - 
      2*lp*X1c[j]*X2c[j]*Y1s[j]*curvature[j] - lp*X1c[j]*X1c[j]*Y2s[j]*curvature[j] + 
      2*lp*Y1c[j]*Y1s[j]*Z20[j]*torsion[j] - lp*X1c[j]*X1c[j]*Z2s[j]*torsion[j] - 
      lp*Y1c[j]*Y1c[j]*Z2s[j]*torsion[j] - lp*Y1s[j]*Y1s[j]*Z2s[j]*torsion[j] - 
      Y1s[j]*Z20[j]*d_X1c_d_varphi[j] - Y1s[j]*Z2c[j]*d_X1c_d_varphi[j] + 
      Y1c[j]*Z2s[j]*d_X1c_d_varphi[j] - X1c[j]*Z2s[j]*d_Y1c_d_varphi[j] + 
      X1c[j]*Z20[j]*d_Y1s_d_varphi[j] - X1c[j]*Z2c[j]*d_Y1s_d_varphi[j] + 
      X1c[j]*Y1s[j]*d_Z2c_d_varphi[j];

    g2s = 4*lp*X2c[j]*Y20[j] - 4*lp*X20[j]*Y2c[j] + 
      lp*X1c[j]*X20[j]*Y1c[j]*curvature[j] - lp*X1c[j]*X2c[j]*Y1c[j]*curvature[j] - 
      2*lp*X1c[j]*X2s[j]*Y1s[j]*curvature[j] - lp*X1c[j]*X1c[j]*Y20[j]*curvature[j] + 
      lp*X1c[j]*X1c[j]*Y2c[j]*curvature[j] - lp*X1c[j]*X1c[j]*Z20[j]*torsion[j] - 
      lp*Y1c[j]*Y1c[j]*Z20[j]*torsion[j] + lp*Y1s[j]*Y1s[j]*Z20[j]*torsion[j] + 
      lp*X1c[j]*X1c[j]*Z2c[j]*torsion[j] + lp*Y1c[j]*Y1c[j]*Z2c[j]*torsion[j] + 
      lp*Y1s[j]*Y1s[j]*Z2c[j]*torsion[j] + Y1c[j]*Z20[j]*d_X1c_d_varphi[j] - 
      Y1c[j]*Z2c[j]*d_X1c_d_varphi[j] - Y1s[j]*Z2s[j]*d_X1c_d_varphi[j] - 
      X1c[j]*Z20[j]*d_Y1c_d_varphi[j] + X1c[j]*Z2c[j]*d_Y1c_d_varphi[j] - 
      X1c[j]*Z2s[j]*d_Y1s_d_varphi[j] + X1c[j]*Y1s[j]*d_Z2s_d_varphi[j];

    if (false) {
      g3s1 = lp*(2*X20[j]*X20[j]*Y1c[j]*curvature[j] + X2c[j]*X2c[j]*Y1c[j]*curvature[j] + X2s[j]*X2s[j]*Y1c[j]*curvature[j] - X1c[j]*X2s[j]*Y2s[j]*curvature[j] + 
		 2*Y1c[j]*Z20[j]*Z20[j]*curvature[j] - 3*Y1c[j]*Z20[j]*Z2c[j]*curvature[j] + Y1c[j]*Z2c[j]*Z2c[j]*curvature[j] - 3*Y1s[j]*Z20[j]*Z2s[j]*curvature[j] + 
		 Y1c[j]*Z2s[j]*Z2s[j]*curvature[j] - 2*Y1c[j]*Y20[j]*Z20[j]*torsion[j] - Y1c[j]*Y2c[j]*Z20[j]*torsion[j] - Y1s[j]*Y2s[j]*Z20[j]*torsion[j] + 
		 4*Y1c[j]*Y20[j]*Z2c[j]*torsion[j] - Y1c[j]*Y2c[j]*Z2c[j]*torsion[j] + 5*Y1s[j]*Y2s[j]*Z2c[j]*torsion[j] - 
		 X1c[j]*X2s[j]*Z2s[j]*torsion[j] + 4*Y1s[j]*Y20[j]*Z2s[j]*torsion[j] - 5*Y1s[j]*Y2c[j]*Z2s[j]*torsion[j] - 
		 Y1c[j]*Y2s[j]*Z2s[j]*torsion[j] - X1c[j]*X2c[j]*(Y20[j]*curvature[j] + Y2c[j]*curvature[j] + (Z20[j] + Z2c[j])*torsion[j]) - 
		 X20[j]*(3*X2c[j]*Y1c[j]*curvature[j] + 3*X2s[j]*Y1s[j]*curvature[j] + 
			 2*X1c[j]*(Y20[j]*curvature[j] - 2*Y2c[j]*curvature[j] + (Z20[j] - 2*Z2c[j])*torsion[j]))) - 2*Y20[j]*Z2c[j]*d_X1c_d_varphi[j] + 
	2*Y1c[j]*Z20[j]*d_X20_d_varphi[j] - 2*Y1c[j]*Z2c[j]*d_X20_d_varphi[j] - 
	2*Y1s[j]*Z2s[j]*d_X20_d_varphi[j] - Y1c[j]*Z20[j]*d_X2c_d_varphi[j] + Y1c[j]*Z2c[j]*d_X2c_d_varphi[j] + 
	Y1s[j]*Z2s[j]*d_X2c_d_varphi[j] - Y1s[j]*Z20[j]*d_X2s_d_varphi[j] - Y1s[j]*Z2c[j]*d_X2s_d_varphi[j] + 
	Y1c[j]*Z2s[j]*d_X2s_d_varphi[j] - 2*X2c[j]*Z20[j]*d_Y1c_d_varphi[j] + 
	2*X20[j]*Z2c[j]*d_Y1c_d_varphi[j] - 2*X2s[j]*Z20[j]*d_Y1s_d_varphi[j] + 
	4*X2s[j]*Z2c[j]*d_Y1s_d_varphi[j] + 2*X20[j]*Z2s[j]*d_Y1s_d_varphi[j] - 
	4*X2c[j]*Z2s[j]*d_Y1s_d_varphi[j] - 2*X1c[j]*Z20[j]*d_Y20_d_varphi[j] + 
	2*X1c[j]*Z2c[j]*d_Y20_d_varphi[j] + X1c[j]*Z20[j]*d_Y2c_d_varphi[j] - X1c[j]*Z2c[j]*d_Y2c_d_varphi[j] - 
	X1c[j]*Z2s[j]*d_Y2s_d_varphi[j] - 2*X20[j]*Y1c[j]*d_Z20_d_varphi[j] + 
	2*X2c[j]*Y1c[j]*d_Z20_d_varphi[j] + 2*X2s[j]*Y1s[j]*d_Z20_d_varphi[j] + 
	2*X1c[j]*Y20[j]*d_Z20_d_varphi[j] + X20[j]*Y1c[j]*d_Z2c_d_varphi[j] - X2c[j]*Y1c[j]*d_Z2c_d_varphi[j] - 
	X2s[j]*Y1s[j]*d_Z2c_d_varphi[j] - X1c[j]*Y20[j]*d_Z2c_d_varphi[j] + 
	Y2c[j]*(2*Z20[j]*d_X1c_d_varphi[j] + X1c[j]*(-2*d_Z20_d_varphi[j] + d_Z2c_d_varphi[j])) - 
	X2s[j]*Y1c[j]*d_Z2s_d_varphi[j] + X20[j]*Y1s[j]*d_Z2s_d_varphi[j] + X2c[j]*Y1s[j]*d_Z2s_d_varphi[j] + 
	X1c[j]*Y2s[j]*d_Z2s_d_varphi[j];

      g3s3 = lp*(-(X2c[j]*X2c[j]*Y1c[j]*curvature[j]) + X2s[j]*X2s[j]*Y1c[j]*curvature[j] - X1c[j]*X2s[j]*Y2s[j]*curvature[j] + Y1c[j]*Z20[j]*Z2c[j]*curvature[j] - 
		 Y1c[j]*Z2c[j]*Z2c[j]*curvature[j] - Y1s[j]*Z20[j]*Z2s[j]*curvature[j] - 2*Y1s[j]*Z2c[j]*Z2s[j]*curvature[j] + Y1c[j]*Z2s[j]*Z2s[j]*curvature[j] - 
		 3*Y1c[j]*Y2c[j]*Z20[j]*torsion[j] + 3*Y1s[j]*Y2s[j]*Z20[j]*torsion[j] + 2*Y1c[j]*Y20[j]*Z2c[j]*torsion[j] + 
		 Y1c[j]*Y2c[j]*Z2c[j]*torsion[j] + Y1s[j]*Y2s[j]*Z2c[j]*torsion[j] - X1c[j]*X2s[j]*Z2s[j]*torsion[j] - 
		 2*Y1s[j]*Y20[j]*Z2s[j]*torsion[j] + Y1s[j]*Y2c[j]*Z2s[j]*torsion[j] - Y1c[j]*Y2s[j]*Z2s[j]*torsion[j] + 
		 X20[j]*(X2c[j]*Y1c[j]*curvature[j] - X2s[j]*Y1s[j]*curvature[j] + 2*X1c[j]*(Y2c[j]*curvature[j] + Z2c[j]*torsion[j])) + 
		 X2c[j]*(-2*X2s[j]*Y1s[j]*curvature[j] + X1c[j]*(-3*Y20[j]*curvature[j] + Y2c[j]*curvature[j] + (-3*Z20[j] + Z2c[j])*torsion[j]))) - 
	2*Y20[j]*Z2c[j]*d_X1c_d_varphi[j] + Y1c[j]*Z20[j]*d_X2c_d_varphi[j] - Y1c[j]*Z2c[j]*d_X2c_d_varphi[j] - 
	Y1s[j]*Z2s[j]*d_X2c_d_varphi[j] - Y1s[j]*Z20[j]*d_X2s_d_varphi[j] - Y1s[j]*Z2c[j]*d_X2s_d_varphi[j] + 
	Y1c[j]*Z2s[j]*d_X2s_d_varphi[j] - 2*X2c[j]*Z20[j]*d_Y1c_d_varphi[j] + 
	2*X20[j]*Z2c[j]*d_Y1c_d_varphi[j] + 2*X2s[j]*Z20[j]*d_Y1s_d_varphi[j] - 
	2*X20[j]*Z2s[j]*d_Y1s_d_varphi[j] - X1c[j]*Z20[j]*d_Y2c_d_varphi[j] + X1c[j]*Z2c[j]*d_Y2c_d_varphi[j] - 
	X1c[j]*Z2s[j]*d_Y2s_d_varphi[j] - X20[j]*Y1c[j]*d_Z2c_d_varphi[j] + X2c[j]*Y1c[j]*d_Z2c_d_varphi[j] + 
	X2s[j]*Y1s[j]*d_Z2c_d_varphi[j] + X1c[j]*Y20[j]*d_Z2c_d_varphi[j] + 
	Y2c[j]*(2*Z20[j]*d_X1c_d_varphi[j] - X1c[j]*d_Z2c_d_varphi[j]) - X2s[j]*Y1c[j]*d_Z2s_d_varphi[j] + 
	X20[j]*Y1s[j]*d_Z2s_d_varphi[j] + X2c[j]*Y1s[j]*d_Z2s_d_varphi[j] + X1c[j]*Y2s[j]*d_Z2s_d_varphi[j];

      g3c1 = -(lp*(2*X20[j]*X20[j]*Y1s[j]*curvature[j] + X2c[j]*X2c[j]*Y1s[j]*curvature[j] + X2s[j]*X2s[j]*Y1s[j]*curvature[j] - X1c[j]*X2s[j]*Y20[j]*curvature[j] - 
		   5*X1c[j]*X2s[j]*Y2c[j]*curvature[j] + 2*Y1s[j]*Z20[j]*Z20[j]*curvature[j] + 3*Y1s[j]*Z20[j]*Z2c[j]*curvature[j] + Y1s[j]*Z2c[j]*Z2c[j]*curvature[j] - 
		   3*Y1c[j]*Z20[j]*Z2s[j]*curvature[j] + Y1s[j]*Z2s[j]*Z2s[j]*curvature[j] - X1c[j]*X2s[j]*Z20[j]*torsion[j] - 
		   2*Y1s[j]*Y20[j]*Z20[j]*torsion[j] + Y1s[j]*Y2c[j]*Z20[j]*torsion[j] - Y1c[j]*Y2s[j]*Z20[j]*torsion[j] - 
		   5*X1c[j]*X2s[j]*Z2c[j]*torsion[j] - 4*Y1s[j]*Y20[j]*Z2c[j]*torsion[j] - Y1s[j]*Y2c[j]*Z2c[j]*torsion[j] - 
		   5*Y1c[j]*Y2s[j]*Z2c[j]*torsion[j] + 4*Y1c[j]*Y20[j]*Z2s[j]*torsion[j] + 5*Y1c[j]*Y2c[j]*Z2s[j]*torsion[j] - 
		   Y1s[j]*Y2s[j]*Z2s[j]*torsion[j] + 5*X1c[j]*X2c[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j]) + 
		   X20[j]*(-3*X2s[j]*Y1c[j]*curvature[j] + 3*X2c[j]*Y1s[j]*curvature[j] + 4*X1c[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j])))) + 
	2*Y20[j]*Z2s[j]*d_X1c_d_varphi[j] + 4*Y2c[j]*Z2s[j]*d_X1c_d_varphi[j] - 
	2*Y1s[j]*Z20[j]*d_X20_d_varphi[j] - 2*Y1s[j]*Z2c[j]*d_X20_d_varphi[j] + 
	2*Y1c[j]*Z2s[j]*d_X20_d_varphi[j] - Y1s[j]*Z20[j]*d_X2c_d_varphi[j] - Y1s[j]*Z2c[j]*d_X2c_d_varphi[j] + 
	Y1c[j]*Z2s[j]*d_X2c_d_varphi[j] + Y1c[j]*Z20[j]*d_X2s_d_varphi[j] - Y1c[j]*Z2c[j]*d_X2s_d_varphi[j] - 
	Y1s[j]*Z2s[j]*d_X2s_d_varphi[j] + 2*X2s[j]*Z20[j]*d_Y1c_d_varphi[j] + 
	4*X2s[j]*Z2c[j]*d_Y1c_d_varphi[j] - 2*X20[j]*Z2s[j]*d_Y1c_d_varphi[j] - 
	4*X2c[j]*Z2s[j]*d_Y1c_d_varphi[j] - 2*X2c[j]*Z20[j]*d_Y1s_d_varphi[j] + 
	2*X20[j]*Z2c[j]*d_Y1s_d_varphi[j] - 2*X1c[j]*Z2s[j]*d_Y20_d_varphi[j] - 
	X1c[j]*Z2s[j]*d_Y2c_d_varphi[j] - X1c[j]*Z20[j]*d_Y2s_d_varphi[j] + X1c[j]*Z2c[j]*d_Y2s_d_varphi[j] - 
	2*X2s[j]*Y1c[j]*d_Z20_d_varphi[j] + 2*X20[j]*Y1s[j]*d_Z20_d_varphi[j] + 
	2*X2c[j]*Y1s[j]*d_Z20_d_varphi[j] - X2s[j]*Y1c[j]*d_Z2c_d_varphi[j] + X20[j]*Y1s[j]*d_Z2c_d_varphi[j] + 
	X2c[j]*Y1s[j]*d_Z2c_d_varphi[j] + Y2s[j]*
	(-2*Z20[j]*d_X1c_d_varphi[j] - 4*Z2c[j]*d_X1c_d_varphi[j] + 
	 X1c[j]*(2*d_Z20_d_varphi[j] + d_Z2c_d_varphi[j])) - X20[j]*Y1c[j]*d_Z2s_d_varphi[j] + 
	X2c[j]*Y1c[j]*d_Z2s_d_varphi[j] + X2s[j]*Y1s[j]*d_Z2s_d_varphi[j] + X1c[j]*Y20[j]*d_Z2s_d_varphi[j] - 
	X1c[j]*Y2c[j]*d_Z2s_d_varphi[j];

      g3c3 = -(lp*(X2c[j]*X2c[j]*Y1s[j]*curvature[j] - X2s[j]*X2s[j]*Y1s[j]*curvature[j] - 3*X1c[j]*X2s[j]*Y20[j]*curvature[j] + X1c[j]*X2s[j]*Y2c[j]*curvature[j] + 
		   Y1s[j]*Z20[j]*Z2c[j]*curvature[j] + Y1s[j]*Z2c[j]*Z2c[j]*curvature[j] + Y1c[j]*Z20[j]*Z2s[j]*curvature[j] - 2*Y1c[j]*Z2c[j]*Z2s[j]*curvature[j] - 
		   Y1s[j]*Z2s[j]*Z2s[j]*curvature[j] - 3*X1c[j]*X2s[j]*Z20[j]*torsion[j] - 3*Y1s[j]*Y2c[j]*Z20[j]*torsion[j] - 
		   3*Y1c[j]*Y2s[j]*Z20[j]*torsion[j] + X1c[j]*X2s[j]*Z2c[j]*torsion[j] + 2*Y1s[j]*Y20[j]*Z2c[j]*torsion[j] - 
		   Y1s[j]*Y2c[j]*Z2c[j]*torsion[j] + Y1c[j]*Y2s[j]*Z2c[j]*torsion[j] + 2*Y1c[j]*Y20[j]*Z2s[j]*torsion[j] + 
		   Y1c[j]*Y2c[j]*Z2s[j]*torsion[j] + Y1s[j]*Y2s[j]*Z2s[j]*torsion[j] + 
		   X2c[j]*(-2*X2s[j]*Y1c[j]*curvature[j] + X1c[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j])) + 
		   X20[j]*(X2s[j]*Y1c[j]*curvature[j] + X2c[j]*Y1s[j]*curvature[j] + 2*X1c[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j])))) + 
	2*Y20[j]*Z2s[j]*d_X1c_d_varphi[j] - Y1s[j]*Z20[j]*d_X2c_d_varphi[j] - Y1s[j]*Z2c[j]*d_X2c_d_varphi[j] + 
	Y1c[j]*Z2s[j]*d_X2c_d_varphi[j] - Y1c[j]*Z20[j]*d_X2s_d_varphi[j] + Y1c[j]*Z2c[j]*d_X2s_d_varphi[j] + 
	Y1s[j]*Z2s[j]*d_X2s_d_varphi[j] + 2*X2s[j]*Z20[j]*d_Y1c_d_varphi[j] - 
	2*X20[j]*Z2s[j]*d_Y1c_d_varphi[j] + 2*X2c[j]*Z20[j]*d_Y1s_d_varphi[j] - 
	2*X20[j]*Z2c[j]*d_Y1s_d_varphi[j] - X1c[j]*Z2s[j]*d_Y2c_d_varphi[j] + X1c[j]*Z20[j]*d_Y2s_d_varphi[j] - 
	X1c[j]*Z2c[j]*d_Y2s_d_varphi[j] - X2s[j]*Y1c[j]*d_Z2c_d_varphi[j] + X20[j]*Y1s[j]*d_Z2c_d_varphi[j] + 
	X2c[j]*Y1s[j]*d_Z2c_d_varphi[j] + Y2s[j]*(-2*Z20[j]*d_X1c_d_varphi[j] + X1c[j]*d_Z2c_d_varphi[j]) + 
	X20[j]*Y1c[j]*d_Z2s_d_varphi[j] - X2c[j]*Y1c[j]*d_Z2s_d_varphi[j] - X2s[j]*Y1s[j]*d_Z2s_d_varphi[j] - 
	X1c[j]*Y20[j]*d_Z2s_d_varphi[j] + X1c[j]*Y2c[j]*d_Z2s_d_varphi[j];

      g40 = -2*(-3*lp*(-((Y2s[j]*Z2c[j] - Y2c[j]*Z2s[j])*(Z20[j]*curvature[j] - Y20[j]*torsion[j])) + 
		       X20[j]*(X2s[j]*(Y2c[j]*curvature[j] + Z2c[j]*torsion[j]) - X2c[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j]))) - 
		2*Y2c[j]*Z2s[j]*d_X20_d_varphi[j] - Y20[j]*Z2s[j]*d_X2c_d_varphi[j] - 
		Y2c[j]*Z20[j]*d_X2s_d_varphi[j] + Y20[j]*Z2c[j]*d_X2s_d_varphi[j] - 
		2*X2s[j]*Z2c[j]*d_Y20_d_varphi[j] + 2*X2c[j]*Z2s[j]*d_Y20_d_varphi[j] - 
		X2s[j]*Z20[j]*d_Y2c_d_varphi[j] + X20[j]*Z2s[j]*d_Y2c_d_varphi[j] + X2c[j]*Z20[j]*d_Y2s_d_varphi[j] - 
		X20[j]*Z2c[j]*d_Y2s_d_varphi[j] + 2*X2s[j]*Y2c[j]*d_Z20_d_varphi[j] + 
		X2s[j]*Y20[j]*d_Z2c_d_varphi[j] + Y2s[j]*
		(2*Z2c[j]*d_X20_d_varphi[j] + Z20[j]*d_X2c_d_varphi[j] - 2*X2c[j]*d_Z20_d_varphi[j] - 
		 X20[j]*d_Z2c_d_varphi[j]) - X2c[j]*Y20[j]*d_Z2s_d_varphi[j] + X20[j]*Y2c[j]*d_Z2s_d_varphi[j]);

      g4s2 = 4*(lp*(Y2c[j]*Z20[j]*Z20[j]*curvature[j] - Y20[j]*Z20[j]*Z2c[j]*curvature[j] - Y2s[j]*Z2c[j]*Z2s[j]*curvature[j] + Y2c[j]*Z2s[j]*Z2s[j]*curvature[j] - 
		    Y20[j]*Y2c[j]*Z20[j]*torsion[j] + Y20[j]*Y20[j]*Z2c[j]*torsion[j] + Y2s[j]*Y2s[j]*Z2c[j]*torsion[j] - Y2c[j]*Y2s[j]*Z2s[j]*torsion[j] - 
		    X20[j]*X2c[j]*(Y20[j]*curvature[j] + Z20[j]*torsion[j]) + X20[j]*X20[j]*(Y2c[j]*curvature[j] + Z2c[j]*torsion[j]) + 
		    X2s[j]*X2s[j]*(Y2c[j]*curvature[j] + Z2c[j]*torsion[j]) - X2c[j]*X2s[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j])) - 
		Y20[j]*Z2c[j]*d_X20_d_varphi[j] - Y2s[j]*Z2c[j]*d_X2s_d_varphi[j] - X2c[j]*Z20[j]*d_Y20_d_varphi[j] + 
		X20[j]*Z2c[j]*d_Y20_d_varphi[j] + X2s[j]*Z2c[j]*d_Y2s_d_varphi[j] - X2c[j]*Z2s[j]*d_Y2s_d_varphi[j] + 
		X2c[j]*Y20[j]*d_Z20_d_varphi[j] + X2c[j]*Y2s[j]*d_Z2s_d_varphi[j] + 
		Y2c[j]*(Z20[j]*d_X20_d_varphi[j] + Z2s[j]*d_X2s_d_varphi[j] - X20[j]*d_Z20_d_varphi[j] - 
			X2s[j]*d_Z2s_d_varphi[j]));

      g4s4 = 2*(lp*(Y2c[j]*Z20[j]*Z2c[j]*curvature[j] - Y20[j]*Z2c[j]*Z2c[j]*curvature[j] - Y2s[j]*Z20[j]*Z2s[j]*curvature[j] + Y20[j]*Z2s[j]*Z2s[j]*curvature[j] - 
		    Y2c[j]*Y2c[j]*Z20[j]*torsion[j] + Y2s[j]*Y2s[j]*Z20[j]*torsion[j] + Y20[j]*Y2c[j]*Z2c[j]*torsion[j] - Y20[j]*Y2s[j]*Z2s[j]*torsion[j] - 
		    X2c[j]*X2c[j]*(Y20[j]*curvature[j] + Z20[j]*torsion[j]) + X2s[j]*X2s[j]*(Y20[j]*curvature[j] + Z20[j]*torsion[j]) + 
		    X20[j]*X2c[j]*(Y2c[j]*curvature[j] + Z2c[j]*torsion[j]) - X20[j]*X2s[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j])) - 
		Y20[j]*Z2c[j]*d_X2c_d_varphi[j] - Y2s[j]*Z20[j]*d_X2s_d_varphi[j] + Y20[j]*Z2s[j]*d_X2s_d_varphi[j] - 
		X2c[j]*Z20[j]*d_Y2c_d_varphi[j] + X20[j]*Z2c[j]*d_Y2c_d_varphi[j] + X2s[j]*Z20[j]*d_Y2s_d_varphi[j] - 
		X20[j]*Z2s[j]*d_Y2s_d_varphi[j] + X2c[j]*Y20[j]*d_Z2c_d_varphi[j] + 
		Y2c[j]*(Z20[j]*d_X2c_d_varphi[j] - X20[j]*d_Z2c_d_varphi[j]) - X2s[j]*Y20[j]*d_Z2s_d_varphi[j] + 
		X20[j]*Y2s[j]*d_Z2s_d_varphi[j]);

      g4c2 = -4*(lp*(Y2s[j]*Z20[j]*Z20[j]*curvature[j] + Y2s[j]*Z2c[j]*Z2c[j]*curvature[j] - Y20[j]*Z20[j]*Z2s[j]*curvature[j] - Y2c[j]*Z2c[j]*Z2s[j]*curvature[j] - 
		     Y20[j]*Y2s[j]*Z20[j]*torsion[j] - Y2c[j]*Y2s[j]*Z2c[j]*torsion[j] + Y20[j]*Y20[j]*Z2s[j]*torsion[j] + Y2c[j]*Y2c[j]*Z2s[j]*torsion[j] - 
		     X20[j]*X2s[j]*(Y20[j]*curvature[j] + Z20[j]*torsion[j]) - X2c[j]*X2s[j]*(Y2c[j]*curvature[j] + Z2c[j]*torsion[j]) + 
		     X20[j]*X20[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j]) + X2c[j]*X2c[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j])) - 
		 Y20[j]*Z2s[j]*d_X20_d_varphi[j] - Y2c[j]*Z2s[j]*d_X2c_d_varphi[j] - X2s[j]*Z20[j]*d_Y20_d_varphi[j] + 
		 X20[j]*Z2s[j]*d_Y20_d_varphi[j] - X2s[j]*Z2c[j]*d_Y2c_d_varphi[j] + X2c[j]*Z2s[j]*d_Y2c_d_varphi[j] + 
		 X2s[j]*Y20[j]*d_Z20_d_varphi[j] + X2s[j]*Y2c[j]*d_Z2c_d_varphi[j] + 
		 Y2s[j]*(Z20[j]*d_X20_d_varphi[j] + Z2c[j]*d_X2c_d_varphi[j] - X20[j]*d_Z20_d_varphi[j] - 
			 X2c[j]*d_Z2c_d_varphi[j]));

      g4c4 = -2*(lp*(Y2s[j]*Z20[j]*Z2c[j]*curvature[j] + Y2c[j]*Z20[j]*Z2s[j]*curvature[j] - 2*Y20[j]*Z2c[j]*Z2s[j]*curvature[j] - 
		     2*Y2c[j]*Y2s[j]*Z20[j]*torsion[j] + Y20[j]*Y2s[j]*Z2c[j]*torsion[j] + Y20[j]*Y2c[j]*Z2s[j]*torsion[j] + 
		     X20[j]*X2s[j]*(Y2c[j]*curvature[j] + Z2c[j]*torsion[j]) + 
		     X2c[j]*(-2*X2s[j]*(Y20[j]*curvature[j] + Z20[j]*torsion[j]) + X20[j]*(Y2s[j]*curvature[j] + Z2s[j]*torsion[j]))) - 
		 Y20[j]*Z2s[j]*d_X2c_d_varphi[j] + Y2c[j]*Z20[j]*d_X2s_d_varphi[j] - Y20[j]*Z2c[j]*d_X2s_d_varphi[j] - 
		 X2s[j]*Z20[j]*d_Y2c_d_varphi[j] + X20[j]*Z2s[j]*d_Y2c_d_varphi[j] - X2c[j]*Z20[j]*d_Y2s_d_varphi[j] + 
		 X20[j]*Z2c[j]*d_Y2s_d_varphi[j] + X2s[j]*Y20[j]*d_Z2c_d_varphi[j] + 
		 Y2s[j]*(Z20[j]*d_X2c_d_varphi[j] - X20[j]*d_Z2c_d_varphi[j]) + X2c[j]*Y20[j]*d_Z2s_d_varphi[j] - 
		 X20[j]*Y2c[j]*d_Z2s_d_varphi[j]);
    }

    // We consider the system sqrt(g) = 0 and
    // d (sqrtg) / d theta = 0.
    // We algebraically eliminate r in "20200322-02 Max r for Garren Boozer.nb", in the section
    // "Keeping first 3 orders in the Jacobian".
    // We end up with the form in "20200322-01 Max r for GarrenBoozer.docx":
    // K0 + K2s*sin(2*theta) + K2c*cos(2*theta) + K4s*sin(4*theta) + K4c*cos(4*theta) = 0.

    K0 = 2*g1c*g1c*g20 - 3*g1c*g1c*g2c + 8*g0*g2c*g2c + 8*g0*g2s*g2s;

    K2s = 2*g1c*g1c*g2s;

    K2c = -2*g1c*g1c*g20 + 2*g1c*g1c*g2c;

    K4s = g1c*g1c*g2s - 16*g0*g2c*g2s;

    K4c = g1c*g1c*g2c - 8*g0*g2c*g2c + 8*g0*g2s*g2s;

    coefficients[0] = 4*(K4c*K4c + K4s*K4s);

    coefficients[1] = 4*(K4s*K2c - K2s*K4c);

    coefficients[2] = K2s*K2s + K2c*K2c - 4*K0*K4c - 4*K4c*K4c - 4*K4s*K4s;

    coefficients[3] = 2*K0*K2s + 2*K4c*K2s - 4*K4s*K2c;

    coefficients[4] = (K0 + K4c)*(K0 + K4c) - K2c*K2c;

    quartic_roots(coefficients, real_parts, imag_parts);

    if (verbose > 1) {
      std::cout << "g0: " << g0 << "  g1c: " << g1c << std::endl;
      std::cout << "g20: " << g20 << "  g2s: " << g2s << "  g2c: " << g2c << std::endl;
      std::cout << "K0: " << K0 << "  K2s: " << K2s << "  K2c: " << K2c << std::endl;
      std::cout << "K4s: " << K4s << "  K4c: " << K4c << std::endl;
      std::cout << "coefficients: "
		<< coefficients[0] << " "
		<< coefficients[1] << " "
		<< coefficients[2] << " "
		<< coefficients[3] << " "
		<< coefficients[4] << std::endl;
      std::cout << "real parts: "
		<< real_parts[0] << " "
		<< real_parts[1] << " "
		<< real_parts[2] << " "
		<< real_parts[3] << std::endl;
      std::cout << "imag parts: "
		<< imag_parts[0] << " "
		<< imag_parts[1] << " "
		<< imag_parts[2] << " "
		<< imag_parts[3] << std::endl;
    }
    
    // Set a default value for rc that is huge to indicate a true solution has not yet been found.
    rc = 1.0e+30;

    for (jr = 0; jr < 4; jr++) { // Loop over the roots of the equation for w.
      // If root is not purely real, skip it.
      if (std::abs(imag_parts[jr]) > 1e-7) {
	if (verbose > 1) std::cout << "Skipping root with jr=" << jr <<
			   " since imag part is" << imag_parts[jr] << std::endl;
	continue;
      }

      sin2theta = real_parts[jr];

      // Discard any roots that have magnitude larger than 1. (I'm not
      // sure this ever happens, but check to be sure.)
      if (std::abs(sin2theta) > 1) {
	if (verbose > 1) std::cout << "Skipping root with jr=" << jr <<
			   " since sin2theta=" << sin2theta << std::endl;
	continue;
      }

      // Determine varpi by checking which choice gives the smaller residual in the K equation
      abs_cos2theta = sqrt(1 - sin2theta * sin2theta);
      residual_if_varpi_plus  = std::abs(K0 + K2s * sin2theta + K2c *   abs_cos2theta 
				    + K4s * 2 * sin2theta *   abs_cos2theta
				    + K4c * (1 - 2 * sin2theta * sin2theta));
      residual_if_varpi_minus = std::abs(K0 + K2s * sin2theta + K2c * (-abs_cos2theta) 
				    + K4s * 2 * sin2theta * (-abs_cos2theta)
				    + K4c * (1 - 2 * sin2theta * sin2theta));

      if (residual_if_varpi_plus > residual_if_varpi_minus) {
	varpi = -1;
      } else {
	varpi = 1;
      }
      cos2theta = varpi * abs_cos2theta;

      // The next few lines give an older method for computing varpi, which has problems in edge cases
      // where w (the root of the quartic polynomial) is very close to +1 or -1, giving varpi
      // not very close to +1 or -1 due to bad loss of precision.
      //
      //varpi_denominator = ((K4s*2*sin2theta + K2c) * sqrt(1 - sin2theta*sin2theta))
      //if (abs(varpi_denominator) < 1e-8) print *,"WARNING////// varpi_denominator=",varpi_denominator
      //varpi = -(K0 + K2s * sin2theta + K4c*(1 - 2*sin2theta*sin2theta)) / varpi_denominator
      //if (abs(varpi*varpi-1) > 1e-3) print *,"WARNING////// abs(varpi*varpi-1) =",abs(varpi*varpi-1)
      //varpi = nint(varpi) // Ensure varpi is exactly either +1 or -1.
      //cos2theta = varpi * sqrt(1 - sin2theta*sin2theta)

      if (verbose > 1) std::cout << "  jr=" << jr << "  sin2theta=" << sin2theta
				 << "  cos2theta=" << cos2theta << std::endl;

      // To get (sin theta, cos theta) from (sin 2 theta, cos 2 theta), we consider two cases to
      // avoid precision loss when cos2theta is added to or subtracted from 1:
      get_cos_from_cos2 = cos2theta > 0;
      if (get_cos_from_cos2) {
	abs_costheta = sqrt(0.5 * (1 + cos2theta));
      } else {
	abs_sintheta = sqrt(0.5 * (1 - cos2theta));
      }
      for(varsigma = -1; varsigma <= 1; varsigma += 2) { // so varsigma will be either -1 or +1.
	if (get_cos_from_cos2) {
	  costheta = varsigma * abs_costheta;
	  sintheta = sin2theta / (2 * costheta);
	} else {
	  sintheta = varsigma * abs_sintheta;
	  costheta = sin2theta / (2 * sintheta);
	}
	if (verbose > 1) std::cout << "    varsigma=" << varsigma << "  costheta=" << costheta
				   << "  sintheta=" << sintheta << " get_cos_from_cos2=" << get_cos_from_cos2
				   << " abs(costheta*costheta + sintheta*sintheta - 1):"
				   << std::abs(costheta*costheta + sintheta*sintheta - 1) << std::endl;

	// Sanity test
	if (std::abs(costheta*costheta + sintheta*sintheta - 1) > sin2_cos2_1_tol) {
	  std::cout << "Error: sintheta=" << sintheta << "  costheta=" << costheta << std::endl;
	  std::cout << "j=" << j << "  jr=" << jr << "  sin2theta=" << sin2theta << "  cos2theta=" << cos2theta << std::endl;
	  std::cout << "abs(costheta*costheta + sintheta*sintheta - 1):" << std::abs(costheta*costheta + sintheta*sintheta - 1)
		    << std::endl;
	  //if (trim(general_option)==general_option_single) stop
	  throw std::runtime_error("sin^2 + cos^2 is far from 1.");
	}

	quadratic_A = g20 + g2s * sin2theta + g2c * cos2theta;
	quadratic_B = costheta * g1c;
	quadratic_C = g0;
	radical = sqrt(quadratic_B * quadratic_B - 4 * quadratic_A * quadratic_C);
	// sign_quadratic = -1 or +1:
	for (sign_quadratic = -1; sign_quadratic <= 1; sign_quadratic += 2) {
	  rr = (-quadratic_B + sign_quadratic * radical) / (2 * quadratic_A); // This is the quadratic formula.
	  residual = -g1c*sintheta + 2*rr*(g2s*cos2theta - g2c*sin2theta); // Residual in the equation d sqrt(g) / d theta = 0.
	  if (verbose > 1) std::cout << "    Quadratic method: rr=" << rr
				     << "  residual=" << residual << std::endl;
	  if ((rr>0 && std::abs(residual) < acceptable_residual)) {
	    if (rr < rc) {// If this is a new minimum
	      rc = rr;
	      sintheta_at_rc = sintheta;
	      costheta_at_rc = costheta;
	      if (verbose > 1) std::cout << "      New minimum: rc =" << rc << std::endl;
	    }
	  }
	}
      } // loop over 2 signs of varsigma
    } // loop over the 4 roots of w polynomial
    r_hat_singularity_robust[j] = rc;
  } // loop over nphi
  
  r_singularity_robust = r_hat_singularity_robust.min();
  
  if (verbose > 0) {
    end_time = std::clock();
    auto end = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for r_singularity from chrono:           "
              << elapsed.count() << std::endl;
    std::cout << "Time for r_singularity from ctime (CPU time): "
              << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
  }
}


    /*
           // Try to get r using the simpler method, the equation that is linear in r.
           denominator = 2*(g2s*cos2theta - g2c*sin2theta)
           //if (abs(denominator) > 1e-8) then // This method cannot be used if we would need to divide by 0
           if (.false.) then
              rr = g1c*sintheta / denominator
              residual = g0 + rr*g1c*costheta + rr*rr*(g20 + g2s*sin2theta + g2c*cos2theta) // Residual in the equation sqrt(g)=0.
              if (local_verbose) print *,"    Linear method: rr=",rr,"  residual=",residual
              if ((rr>0) .and. (abs(residual) < 1e-5)) then
                 if (rr < rc) then // If this is a new minimum
                    rc = rr
                    sintheta_at_rc = sintheta
                    costheta_at_rc = costheta
                    if (local_verbose) print *,"      New minimum: rc =",rc
                 end if
              end if
           else
              // Use the more complicated method to determine rr by solving a quadratic equation.
              quadratic_A = g20 + g2s * sin2theta + g2c * cos2theta
              quadratic_B = costheta * g1c
              quadratic_C = g0
              radical = sqrt(quadratic_B * quadratic_B - 4 * quadratic_A * quadratic_C)
              do sign_quadratic = -1, 1, 2 // So sign_quadratic = +1 or -1
                 rr = (-quadratic_B + sign_quadratic * radical) / (2 * quadratic_A) // This is the quadratic formula.
                 residual = -g1c*sintheta + 2*rr*(g2s*cos2theta - g2c*sin2theta) // Residual in the equation d sqrt(g) / d theta = 0.
                 if (local_verbose) print *,"    Quadratic method: rr=",rr,"  residual=",residual
                 if ((rr>0) .and. (abs(residual) < 1e-5)) then
                    if (rr < rc) then // If this is a new minimum
                       rc = rr
                       sintheta_at_rc = sintheta
                       costheta_at_rc = costheta
                       if (local_verbose) print *,"      New minimum: rc =",rc
                    end if
                 end if
              end do
           end if
        end do
     end do
     r_singularity_basic_vs_zeta(j) = rc
    */
