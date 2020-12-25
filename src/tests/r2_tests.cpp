#include <vector>
#include <string>
#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

/** Verify that the O(r^2) equations have been solved, using a
    different method than the method used to solve them originally.
 */
TEST_CASE("Check O(r^2) equations") {
  Vector fX0, fXs, fXc, fY0, fYs, fYc, residual1, residual2;
  std::vector<std::string> configs = {
    "r2 section 5.1",
    "r2 section 5.2",
    "r2 section 5.3",
    "r2 section 5.4",
    "r2 section 5.5"};

  qscfloat half = 0.5, tol;
  
  for (int jconfig = 0; jconfig < configs.size(); jconfig++) {
    CAPTURE(jconfig);
    Qsc q(configs[jconfig]);
    if (single) {
      tol = 2.0e-4;
    } else {
      tol = 1.0e-12;
    }
    
    fX0.resize(q.nphi, 0.0);
    fXs.resize(q.nphi, 0.0);
    fXc.resize(q.nphi, 0.0);
    fY0.resize(q.nphi, 0.0);
    fYs.resize(q.nphi, 0.0);
    fYc.resize(q.nphi, 0.0);
    
    fX0 = q.d_X20_d_varphi - q.torsion * q.abs_G0_over_B0 * q.Y20 + q.curvature * q.abs_G0_over_B0 * q.Z20
      -4*q.sG*q.spsi*q.abs_G0_over_B0*(q.Y2c * q.Z2s - q.Y2s * q.Z2c)
      - q.spsi * q.I2_over_B0 * (half * q.curvature * q.X1c * q.Y1c - 2 * q.Y20) * q.abs_G0_over_B0 + q.abs_G0_over_B0 * q.beta_1s * q.Y1c * half;
    
    fXs = q.d_X2s_d_varphi - 2 * q.iota_N * q.X2c - q.torsion * q.abs_G0_over_B0 * q.Y2s + q.curvature * q.abs_G0_over_B0 * q.Z2s
      -4*q.sG*q.spsi*q.abs_G0_over_B0*(-q.Y20 * q.Z2c + q.Y2c * q.Z20)
      -q.spsi * q.I2_over_B0 * (half * q.curvature * q.X1c * q.Y1s - 2 * q.Y2s) * q.abs_G0_over_B0 - q.abs_G0_over_B0 * q.beta_1s * q.Y1s * half;
    
    fXc = q.d_X2c_d_varphi + 2 * q.iota_N * q.X2s - q.torsion * q.abs_G0_over_B0 * q.Y2c + q.curvature * q.abs_G0_over_B0 * q.Z2c
      -4*q.sG*q.spsi*q.abs_G0_over_B0*(q.Y20 * q.Z2s - q.Y2s * q.Z20)
      -q.spsi * q.I2_over_B0 * (half * q.curvature * q.X1c * q.Y1c - 2 * q.Y2c) * q.abs_G0_over_B0 - q.abs_G0_over_B0 * q.beta_1s * q.Y1c * half;
    
    fY0 = q.d_Y20_d_varphi + q.torsion * q.abs_G0_over_B0 * q.X20 - 4*q.sG*q.spsi*q.abs_G0_over_B0*(q.X2s * q.Z2c - q.X2c * q.Z2s)
      -q.spsi * q.I2_over_B0 * (-half * q.curvature*q.X1c*q.X1c + 2*q.X20) * q.abs_G0_over_B0 - q.abs_G0_over_B0 * q.beta_1s * q.X1c * half;
    
    fYs = q.d_Y2s_d_varphi - 2 * q.iota_N * q.Y2c + q.torsion * q.abs_G0_over_B0 * q.X2s
      -4*q.sG*q.spsi*q.abs_G0_over_B0*(q.X20 * q.Z2c - q.X2c * q.Z20) - 2*q.spsi* q.I2_over_B0 * q.X2s * q.abs_G0_over_B0;
    
    fYc = q.d_Y2c_d_varphi + 2 * q.iota_N * q.Y2s + q.torsion * q.abs_G0_over_B0 * q.X2c
      -4*q.sG*q.spsi*q.abs_G0_over_B0*(q.X2s * q.Z20 - q.X20 * q.Z2s)
      -q.spsi * q.I2_over_B0 * (-half * q.curvature * q.X1c * q.X1c + 2 * q.X2c) * q.abs_G0_over_B0 + q.abs_G0_over_B0 * q.beta_1s * q.X1c * half;
    
    // Check equations (A41)-(A42) in Landreman & Sengupta, J Plasma Physics (2019):
    residual1 = q.X1c * fXs - q.Y1s * fY0 + q.Y1c * fYs - q.Y1s * fYc;
    residual2 = -q.X1c * fX0 + q.X1c * fXc - q.Y1c * fY0 + q.Y1s * fYs + q.Y1c * fYc;
    for (int j = 0; j < q.nphi; j++) {
      CHECK(Approx(residual1[j]).epsilon(tol) == 0.0);
      CHECK(Approx(residual2[j]).epsilon(tol) == 0.0);
    }
    std::cout << "Residual1: " << residual1 << std::endl;
    std::cout << "Residual2: " << residual1 << std::endl;
    
    // Now check equations (A32)-(A33) in Landreman & Sengupta, J Plasma Physics (2019):
    residual1 = -q.X1c * q.Y2c + q.X1c * q.Y20 + q.X2s * q.Y1s + q.X2c * q.Y1c - q.X20 * q.Y1c;
    residual2 = q.X1c * q.Y2s + q.X2c * q.Y1s - q.X2s * q.Y1c + q.X20 * q.Y1s + q.sG * q.spsi * q.X1c * q.curvature * half;
    for (int j = 0; j < q.nphi; j++) {
      CHECK(Approx(residual1[j]).epsilon(tol) == 0.0);
      CHECK(Approx(residual2[j]).epsilon(tol) == 0.0);
    }
    std::cout << "Residual1: " << residual1 << std::endl;
    std::cout << "Residual2: " << residual1 << std::endl;
  }
}

