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

  for (j = 0; j < nphi; j++) {
     // Write sqrt(g) = r * [g0 + r*g1c*cos(theta) + (r^2)*(g20 + g2s*sin(2*theta) + g2c*cos(2*theta) + ...]
     // The coefficients are evaluated in "20200322-02 Max r for Garren Boozer.nb", in the section "Order r^2 construction, quasisymmetry"

    g0 = lp * X1c[j] * Y1s[j];

     //g1s = -2*X20[j]*Y1c[j] + 2*X2c[j]*Y1c[j] + 2*X2s[j]*Y1s[j] + 2*X1c[j]*Y20[j] - 2*X1c[j]*Y2c[j]
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
  }
}

