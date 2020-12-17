#include <iostream>
#include "qsc.hpp"

using namespace qsc;

/** Initialize the axis shape, curvature, and torsion.
 */
void Qsc::init_axis() {
  int j, n;

  // Initialize the phi grid.
  d_phi = 2 * pi / (nfp * nphi);
  phi[0] = 0.0;
  for (j = 1; j < nphi; j++) phi[j] = phi[j - 1] + d_phi;

  // Initialize the axis shape.
  R0 = R0c[0];
  Z0 = Z0c[0];
  R0p = 0.0;
  Z0p = 0.0;
  R0pp = 0.0;
  Z0pp = 0.0;
  R0ppp = 0.0;
  Z0ppp = 0.0;
  for (n = 1; n <= R0c.size(); n++) {
    sinangle = sin(double(n * nfp) * phi);
    cosangle = cos(double(n * nfp) * phi);
    R0 = R0 + R0c[n] * cosangle + R0s[n] * sinangle;
    Z0 = Z0 + Z0c[n] * cosangle + Z0s[n] * sinangle;
    //R0_extended = R0_extended + R0c[n] * cosangle + R0s[n] * sinangle;
    //Z0_extended = Z0_extended + R0c[n] * cosangle + R0s[n] * sinangle;
    R0p = R0p + R0c[n] * (-n*nfp)*sinangle + R0s[n] * (n*nfp)*cosangle;
    Z0p = Z0p + Z0c[n] * (-n*nfp)*sinangle + Z0s[n] * (n*nfp)*cosangle;
    R0pp = R0pp + R0c[n] * (-n*nfp*n*nfp)*cosangle
      + R0s[n] * (-n*nfp*n*nfp)*sinangle;
    Z0pp = Z0pp + Z0c[n] * (-n*nfp*n*nfp)*cosangle
      + Z0s[n] * (-n*nfp*n*nfp)*sinangle;
    R0ppp = R0ppp + R0c[n] * (n*nfp*n*nfp*n*nfp)*sinangle
      + R0s[n] * (-n*nfp*n*nfp*n*nfp)*cosangle;
    Z0ppp = Z0ppp + Z0c[n] * (n*nfp*n*nfp*n*nfp)*sinangle
      + Z0s[n] * (-n*nfp*n*nfp*n*nfp)*cosangle;
  }
  d_l_d_phi = sqrt(R0 * R0 + R0p * R0p + Z0p * Z0p);
  d2_l_d_phi2 = (R0 * R0p + R0p * R0pp + Z0p * Z0pp) / d_l_d_phi;
  B0_over_abs_G0 = nphi / d_l_d_phi.sum();
  abs_G0_over_B0 = 1 / B0_over_abs_G0;
  d_l_d_varphi = abs_G0_over_B0;
  G0 = sG * abs_G0_over_B0 * B0;
  
  std::cout << "R0c:" << R0c << std::endl;
  std::cout << "R0s:" << R0s << std::endl;
  std::cout << "Z0c:" << Z0c << std::endl;
  std::cout << "Z0s:" << Z0s << std::endl;
  std::cout << "phi:" << phi << std::endl;
  std::cout << "R0:" << R0 << std::endl;
}
