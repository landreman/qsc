#include "qsc.hpp"

using namespace qsc;

/** Initialize the axis shape, curvature, and torsion.
 */
void Qsc::init_axis() {
  int j, n;

  // Initialize the phi grid.
  d_phi = 2 * pi / nfp;
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
    sinangle = sin(phi);
    cosangle = cos(phi);
    R0 = R0 + R0c[n] * cosangle + R0s[n] * sinangle;
    Z0 = Z0 + Z0c[n] * cosangle + Z0s[n] * sinangle;
    //R0_extended = R0_extended + R0c[n] * cosangle + R0s[n] * sinangle;
    //Z0_extended = Z0_extended + R0c[n] * cosangle + R0s[n] * sinangle;
    R0p = R0p + R0c[n] * (-n*nfp)*sinangle + R0s[n] * (n*nfp)*cosangle;
    Z0p = Z0p + Z0c[n] * (-n*nfp)*sinangle + Z0s[n] * (n*nfp)*cosangle;
    R0pp = R0pp + R0c[n] * (-n*nfp*n*nfp)*cosangle + R0s[n] * (-n*nfp*n*nfp)*sinangle;
    Z0pp = Z0pp + Z0c[n] * (-n*nfp*n*nfp)*cosangle + Z0s[n] * (-n*nfp*n*nfp)*sinangle;
    R0ppp = R0ppp + R0c[n] * (n*nfp*n*nfp*n*nfp)*sinangle + R0s[n] * (-n*nfp*n*nfp*n*nfp)*cosangle;
    Z0ppp = Z0ppp + Z0c[n] * (n*nfp*n*nfp*n*nfp)*sinangle + Z0s[n] * (-n*nfp*n*nfp*n*nfp)*cosangle;

  }
  
}
