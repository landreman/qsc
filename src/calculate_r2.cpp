#include <algorithm>
#include <cmath>
#include <chrono>
#include "qsc.hpp"

using namespace qsc;

/** Solve the O(r^2) equations for X2, Y2, Z2, and B20.
 */
void Qsc::calculate_r2() {
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) start = std::chrono::steady_clock::now();
  
  if (std::abs(iota_N) < 1.0e-8)
    std::cerr <<
      "Warning: |iota_N| is very small so O(r^2) solve will be poorly conditioned. iota_N="
	      << iota_N << std::endl;

  if (verbose > 0) std::cout << "Beginning O(r^2) calculation" << std::endl;

  qscfloat half = 0.5, quarter = 0.25;
  int j, k;

  G2 = -mu0 * p2 * G0 / (B0 * B0) - iota * I2;
    
  V1 = X1c * X1c + Y1c * Y1c + Y1s * Y1s;
  V2 = 2 * Y1s * Y1c;
  V3 = X1c * X1c + Y1c * Y1c - Y1s * Y1s;

  qscfloat factor = - B0_over_abs_G0 / 8.0;
  matrix_vector_product(d_d_varphi, V1, Z20);
  matrix_vector_product(d_d_varphi, V2, Z2s);
  matrix_vector_product(d_d_varphi, V3, Z2c);
  Z20 *= factor;
  Z2s = factor * (Z2s - 2 * iota_N * V3);
  Z2c = factor * (Z2c + 2 * iota_N * V2);

  matrix_vector_product(d_d_varphi, Z20, d_Z20_d_varphi);
  matrix_vector_product(d_d_varphi, Z2s, d_Z2s_d_varphi);
  matrix_vector_product(d_d_varphi, Z2c, d_Z2c_d_varphi);

  qs = -iota_N * X1c - abs_G0_over_B0 * Y1s * torsion;
  
  matrix_vector_product(d_d_varphi, X1c, qc);
  qc -= abs_G0_over_B0 * Y1c * torsion;

  matrix_vector_product(d_d_varphi, Y1s, rs);
  rs -= iota_N * Y1c;

  matrix_vector_product(d_d_varphi, Y1c, rc);
  rc += iota_N * Y1s + X1c * torsion * abs_G0_over_B0;

  matrix_vector_product(d_d_varphi, Z2s, X2s);
  X2s = B0_over_abs_G0 * (X2s - (2 * iota_N) * Z2c + B0_over_abs_G0 * ( abs_G0_over_B0 * abs_G0_over_B0 * B2s / B0 + (qc * qs + rc * rs) * half)) / curvature;

  matrix_vector_product(d_d_varphi, Z2c, X2c);
  X2c = B0_over_abs_G0 * (X2c + (2 * iota_N) * Z2s - B0_over_abs_G0 * (-abs_G0_over_B0 * abs_G0_over_B0 * B2c / B0 + abs_G0_over_B0 * abs_G0_over_B0 * eta_bar * eta_bar / 2.0 - (qc * qc - qs * qs + rc * rc - rs * rs) * 0.25)) / curvature;

  matrix_vector_product(d_d_varphi, X2s, d_X2s_d_varphi);
  matrix_vector_product(d_d_varphi, X2c, d_X2c_d_varphi);
  
  beta_1s = -4 * spsi * sG * mu0 * p2 * eta_bar * abs_G0_over_B0 / (iota_N * B0 * B0);

  Y2s_from_X20 = -sG * spsi * curvature * curvature / (eta_bar * eta_bar);
  Y2s_inhomogeneous = sG * spsi * (-curvature/2 + curvature*curvature/(eta_bar*eta_bar)*(-X2c + X2s * sigma));

  Y2c_from_X20 = -sG * spsi * curvature * curvature * sigma / (eta_bar * eta_bar);
  Y2c_inhomogeneous = sG * spsi * curvature * curvature / (eta_bar * eta_bar) * (X2s + X2c * sigma);

  /* Note: in the fX* and fY* quantities below, I've omitted the
     contributions from X20 and Y20 to the d/dzeta terms. These
     contributions are handled later when we assemble the large
     matrix.
  */
  fX0_from_X20 = -4 * sG * spsi * abs_G0_over_B0 * (Y2c_from_X20 * Z2s - Y2s_from_X20 * Z2c);
  fX0_from_Y20 = -torsion * abs_G0_over_B0 - 4 * sG * spsi * abs_G0_over_B0 * (Z2s)
    - spsi * I2_over_B0 * (-2) * abs_G0_over_B0;
  fX0_inhomogeneous = curvature * abs_G0_over_B0 * Z20 - 4 * sG * spsi * abs_G0_over_B0 * (Y2c_inhomogeneous * Z2s - Y2s_inhomogeneous * Z2c) 
    - (spsi * I2_over_B0 * half * sG * spsi * abs_G0_over_B0) * curvature + beta_1s * abs_G0_over_B0 / 2 * Y1c;

  fXs_from_X20 = -torsion * abs_G0_over_B0 * Y2s_from_X20 - 4 * spsi * sG * abs_G0_over_B0 * (Y2c_from_X20 * Z20) 
    - spsi * I2_over_B0 * (- 2 * Y2s_from_X20) * abs_G0_over_B0;
  fXs_from_Y20 = - 4 * spsi * sG * abs_G0_over_B0 * (-Z2c + Z20);
  fXs_inhomogeneous = d_X2s_d_varphi - 2 * iota_N * X2c - torsion * abs_G0_over_B0 * Y2s_inhomogeneous + curvature * abs_G0_over_B0 * Z2s
    - 4 * spsi * sG * abs_G0_over_B0 * (Y2c_inhomogeneous * Z20)
    - spsi * I2_over_B0 * ((half * spsi * sG) * curvature - 2 * Y2s_inhomogeneous) * abs_G0_over_B0
    - (half) * abs_G0_over_B0 * beta_1s * Y1s;

  fXc_from_X20 = - torsion * abs_G0_over_B0 * Y2c_from_X20 - 4 * spsi * sG * abs_G0_over_B0 * (-Y2s_from_X20 * Z20)
    - spsi * I2_over_B0 * (- 2 * Y2c_from_X20) * abs_G0_over_B0;
  fXc_from_Y20 = - torsion * abs_G0_over_B0 - 4 * spsi * sG * abs_G0_over_B0 * (Z2s)
    - spsi * I2_over_B0 * (-2) * abs_G0_over_B0;
  fXc_inhomogeneous = d_X2c_d_varphi + 2 * iota_N * X2s - torsion * abs_G0_over_B0 * Y2c_inhomogeneous + curvature * abs_G0_over_B0 * Z2c
    - 4 * spsi * sG * abs_G0_over_B0 * (-Y2s_inhomogeneous * Z20)
    - spsi * I2_over_B0 * ((half * sG * spsi) * curvature - 2 * Y2c_inhomogeneous) * abs_G0_over_B0
    - (half) * abs_G0_over_B0 * beta_1s * Y1c;

  fY0_from_X20 = torsion * abs_G0_over_B0 - spsi * I2_over_B0 * (2) * abs_G0_over_B0;
  fY0_from_Y20 = 0; // Could save a little time here?
  fY0_inhomogeneous = -4 * spsi * sG * abs_G0_over_B0 * (X2s * Z2c - X2c * Z2s)
    - spsi * I2_over_B0 * (-half * curvature * X1c * X1c) * abs_G0_over_B0 - (half) * abs_G0_over_B0 * beta_1s * X1c;

  fYs_from_X20 = -2 * iota_N * Y2c_from_X20 - 4 * spsi * sG * abs_G0_over_B0 * (Z2c);
  fYs_from_Y20 = -2 * iota_N; // Note this is independent of phi
  matrix_vector_product(d_d_varphi, Y2s_inhomogeneous, work1);
  fYs_inhomogeneous = work1 - 2 * iota_N * Y2c_inhomogeneous + torsion * abs_G0_over_B0 * X2s
    - 4 * spsi * sG * abs_G0_over_B0 * (-X2c * Z20) - 2 * spsi * I2_over_B0 * X2s * abs_G0_over_B0;

  fYc_from_X20 = 2 * iota_N * Y2s_from_X20 - 4 * spsi * sG * abs_G0_over_B0 * (-Z2s);
  fYc_from_Y20 = 0; // Could save a little time here?
  matrix_vector_product(d_d_varphi, Y2c_inhomogeneous, work1);
  fYc_inhomogeneous = work1 + 2 * iota_N * Y2s_inhomogeneous + torsion * abs_G0_over_B0 * X2c
    - 4 * spsi * sG * abs_G0_over_B0 * (X2s * Z20)
    - spsi * I2_over_B0 * (-half * curvature * X1c * X1c + 2 * X2c) * abs_G0_over_B0 + half * abs_G0_over_B0 * beta_1s * X1c;

  for (j = 0; j < nphi; j++) {
    for (k = 0; k < nphi; k++) {
      // Handle the terms involving d X_0 / d zeta and d Y_0 / d zeta:
      // ----------------------------------------------------------------

      // Equation 1, terms involving X0:
      // Contributions arise from Y1c * fYs - Y1s * fYc.
      r2_matrix(j, k) = Y1c[j] * d_d_varphi(j, k) * Y2s_from_X20[k]
	- Y1s[j] * d_d_varphi(j, k) * Y2c_from_X20[k];

      // Equation 1, terms involving Y0:
      // Contributions arise from -Y1s * fY0 - Y1s * fYc, and they happen to be equal.
      r2_matrix(j, k + nphi) = -2 * Y1s[j] * d_d_varphi(j, k);

      // Equation 2, terms involving X0:
      // Contributions arise from -X1c * fX0 + Y1s * fYs + Y1c * fYc
      r2_matrix(j + nphi, k) = -X1c[j] * d_d_varphi(j, k)
	+ Y1s[j] * d_d_varphi(j, k) * Y2s_from_X20[k]
	+ Y1c[j] * d_d_varphi(j, k) * Y2c_from_X20[k];

      // Equation 2, terms involving Y0:
      // Contributions arise from -Y1c * fY0 + Y1c * fYc, but they happen to cancel.
      r2_matrix(j + nphi, k + nphi) = 0.0;
    }

    // Now handle the terms involving X_0 and Y_0 without d / d varphi derivatives:
    // ----------------------------------------------------------------

    r2_matrix(j, j       ) = r2_matrix(j, j       )
      + X1c[j] * fXs_from_X20[j] - Y1s[j] * fY0_from_X20[j]
      + Y1c[j] * fYs_from_X20[j] - Y1s[j] * fYc_from_X20[j];
    
    r2_matrix(j, j + nphi) = r2_matrix(j, j + nphi)
      + X1c[j] * fXs_from_Y20[j] - Y1s[j] * fY0_from_Y20[j]
      + Y1c[j] * fYs_from_Y20[j] - Y1s[j] * fYc_from_Y20[j];

    r2_matrix(j + nphi, j       ) = r2_matrix(j + nphi, j      )
      - X1c[j] * fX0_from_X20[j] + X1c[j] * fXc_from_X20[j] - Y1c[j] * fY0_from_X20[j]
      + Y1s[j] * fYs_from_X20[j] + Y1c[j] * fYc_from_X20[j];
    
    r2_matrix(j + nphi, j + nphi) = r2_matrix(j + nphi, j + nphi)
      - X1c[j] * fX0_from_Y20[j] + X1c[j] * fXc_from_Y20[j] - Y1c[j] * fY0_from_Y20[j]
      + Y1s[j] * fYs_from_Y20[j] + Y1c[j] * fYc_from_Y20[j];
  }

  // Assemble the right-hand side
  work1 = -(X1c * fXs_inhomogeneous - Y1s * fY0_inhomogeneous + Y1c * fYs_inhomogeneous - Y1s * fYc_inhomogeneous);
  work2 = -(- X1c * fX0_inhomogeneous + X1c * fXc_inhomogeneous - Y1c * fY0_inhomogeneous + Y1s * fYs_inhomogeneous + Y1c * fYc_inhomogeneous);
  for (j = 0; j < nphi; j++) {
    r2_rhs[j] = work1[j];
    r2_rhs[j + nphi] = work2[j];
  }

  std::chrono::time_point<std::chrono::steady_clock> solve_start, solve_end;
  if (verbose > 0) solve_start = std::chrono::steady_clock::now();
  // Here is the main solve:
  linear_solve(r2_matrix, r2_rhs, r2_ipiv);
  if (verbose > 0) solve_end = std::chrono::steady_clock::now();

  // Extract X20 and Y20 from the solution:
  for (j = 0; j < nphi; j++) {
    X20[j] = r2_rhs[j];
    Y20[j] = r2_rhs[j + nphi];
  }

  // Now that we have X20 and Y20 explicitly, we can reconstruct Y2s, Y2c, and B20:
  Y2s = Y2s_inhomogeneous + Y2s_from_X20 * X20;
  Y2c = Y2c_inhomogeneous + Y2c_from_X20 * X20 + Y20;

  matrix_vector_product(d_d_varphi, X20, d_X20_d_varphi);
  matrix_vector_product(d_d_varphi, Y20, d_Y20_d_varphi);
  matrix_vector_product(d_d_varphi, Y2s, d_Y2s_d_varphi);
  matrix_vector_product(d_d_varphi, Y2c, d_Y2c_d_varphi);
  matrix_vector_product(d_d_varphi, d_X1c_d_varphi, d2_X1c_d_varphi2);
  matrix_vector_product(d_d_varphi, d_Y1c_d_varphi, d2_Y1c_d_varphi2);
  matrix_vector_product(d_d_varphi, d_Y1s_d_varphi, d2_Y1s_d_varphi2);
  matrix_vector_product(d_d_varphi, curvature, d_curvature_d_varphi);
  matrix_vector_product(d_d_varphi, torsion, d_torsion_d_varphi);

  B20 = B0 * (curvature * X20 - B0_over_abs_G0 * d_Z20_d_varphi + half * eta_bar * eta_bar - mu0 * p2 / (B0 * B0)
	      - quarter * B0_over_abs_G0 * B0_over_abs_G0 * (qc * qc + qs * qs + rc * rc + rs * rs));

  B20_grid_variation = B20.max() - B20.min();
  
  qscfloat normalizer = 1.0 / d_l_d_phi.sum();
  work1 = B20 * d_l_d_phi;
  B20_mean = work1.sum() * normalizer;
  B20_anomaly = B20 - B20_mean;
  
  work1 = B20_anomaly * B20_anomaly * d_l_d_phi;
  B20_residual = sqrt(work1.sum() * normalizer) / B0;

  work1 = std::abs(X20);
  grid_max_XY2 = work1.max();
  work1 = std::abs(X2s);
  grid_max_XY2 = std::max(grid_max_XY2, work1.max());
  work1 = std::abs(X2c);
  grid_max_XY2 = std::max(grid_max_XY2, work1.max());
  work1 = std::abs(Y20);
  grid_max_XY2 = std::max(grid_max_XY2, work1.max());
  work1 = std::abs(Y2s);
  grid_max_XY2 = std::max(grid_max_XY2, work1.max());
  work1 = std::abs(Y2c);
  grid_max_XY2 = std::max(grid_max_XY2, work1.max());
  
  work1 = std::abs(Z20);
  grid_max_Z2 = work1.max();
  work1 = std::abs(Z2s);
  grid_max_Z2 = std::max(grid_max_Z2, work1.max());
  work1 = std::abs(Z2c);
  grid_max_Z2 = std::max(grid_max_Z2, work1.max());

  work1 = std::abs(d_X20_d_varphi);
  grid_max_d_XY2_d_varphi = work1.max();
  work1 = std::abs(d_X2s_d_varphi);
  grid_max_d_XY2_d_varphi = std::max(grid_max_d_XY2_d_varphi, work1.max());
  work1 = std::abs(d_X2c_d_varphi);
  grid_max_d_XY2_d_varphi = std::max(grid_max_d_XY2_d_varphi, work1.max());
  work1 = std::abs(d_Y20_d_varphi);
  grid_max_d_XY2_d_varphi = std::max(grid_max_d_XY2_d_varphi, work1.max());
  work1 = std::abs(d_Y2s_d_varphi);
  grid_max_d_XY2_d_varphi = std::max(grid_max_d_XY2_d_varphi, work1.max());
  work1 = std::abs(d_Y2c_d_varphi);
  grid_max_d_XY2_d_varphi = std::max(grid_max_d_XY2_d_varphi, work1.max());
  
  work1 = std::abs(d_Z20_d_varphi);
  grid_max_d_Z2_d_varphi = work1.max();
  work1 = std::abs(d_Z2s_d_varphi);
  grid_max_d_Z2_d_varphi = std::max(grid_max_d_Z2_d_varphi, work1.max());
  work1 = std::abs(d_Z2c_d_varphi);
  grid_max_d_Z2_d_varphi = std::max(grid_max_d_Z2_d_varphi, work1.max());
  
  if (order_r2p1) calculate_r2p1();
  
  ////////////////////////////////////////////////////////////
  
  if (verbose > 0) {
    auto end = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    std::chrono::duration<double> solve_elapsed = solve_end - solve_start;
    std::cout << "Time for calculate_r2: "
              << elapsed.count() << " seconds, " << solve_elapsed.count()
	      << " for solve" << std::endl;
  }
}

/////////////////////////////////////////

void Qsc::r2_diagnostics() {
  std::chrono::time_point<std::chrono::steady_clock> diag_start, diag_end;
  if (verbose > 0) diag_start = std::chrono::steady_clock::now();

  mercier();
  calculate_grad_grad_B_tensor();
  calculate_r_singularity();
  
  if (verbose > 0) {
    diag_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diag_elapsed = diag_end - diag_start;
    std::cout << "Time for r2_diagnostics: "
              << diag_elapsed.count() << std::endl;
  }

}

/////////////////////////////////////////

/** In this subroutine, equations (3.12) and (3.14)-(3.15) of
    Landreman & Sengupta (2019) are implemented. Just enough of the
    O(r^3) terms in X and Y are computed in order to make the m=0 mode
    of B accurate to O(r^2).
 */
void Qsc::calculate_r2p1() {
  // Eq (3.15): lambda = - Q / (2 * sign_G * spsi)
  lambda_for_XY3 = - 0.5 / (sG * spsi) * (-spsi * B0 * abs_G0_over_B0 / (2*G0*G0) * (iota_N * I2 + mu0 * p2 * G0 / (B0 * B0)) + 2 * (X2c * Y2s - X2s * Y2c) 
					  + spsi * B0 / (2*G0) * (abs_G0_over_B0 * X20 * curvature - d_Z20_d_varphi) 
					  + I2 / (4 * G0) * (-abs_G0_over_B0 * torsion * (X1c*X1c + Y1s*Y1s + Y1c*Y1c) + Y1c * d_X1c_d_varphi - X1c * d_Y1c_d_varphi));  
  // The expression above is derived in
  // "20190318-01 Wrick's streamlined Garren-Boozer method, MHD.nb" in
  // the section "Not assuming quasisymmetry".  Note Q = (1/2) *
  // (XYEquation0 without X3 and Y3 terms) where XYEquation0 is the
  // quantity in the above notebook.
  
  // Eq (3.14):
  X3c1 = X1c * lambda_for_XY3;
  Y3c1 = Y1c * lambda_for_XY3;
  Y3s1 = Y1s * lambda_for_XY3;

  matrix_vector_product(d_d_varphi, X3c1, d_X3c1_d_varphi);  
  matrix_vector_product(d_d_varphi, Y3c1, d_Y3c1_d_varphi);  
  matrix_vector_product(d_d_varphi, Y3s1, d_Y3s1_d_varphi);
  
  work1 = std::abs(X3c1);
  grid_max_XY3 = work1.max();
  work1 = std::abs(Y3c1);
  grid_max_XY3 = std::max(grid_max_XY3, work1.max());
  work1 = std::abs(Y3s1);
  grid_max_XY3 = std::max(grid_max_XY3, work1.max());

  work1 = std::abs(d_X3c1_d_varphi);
  grid_max_d_XY3_d_varphi = work1.max();
  work1 = std::abs(d_Y3c1_d_varphi);
  grid_max_d_XY3_d_varphi = std::max(grid_max_d_XY3_d_varphi, work1.max());
  work1 = std::abs(d_Y3s1_d_varphi);
  grid_max_d_XY3_d_varphi = std::max(grid_max_d_XY3_d_varphi, work1.max());
}
