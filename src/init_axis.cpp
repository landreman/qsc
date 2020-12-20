#include <ctime>
#include <chrono>
#include <iostream>
#include "qsc.hpp"

using namespace qsc;

/** Initialize the axis shape, curvature, and torsion.
 */
void Qsc::init_axis() {
  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }

  int j, k, n;

  // Initialize the phi grid.
  d_phi = 2 * pi / (nfp * nphi);
  phi[0] = 0.0;
  for (j = 1; j < nphi; j++) phi[j] = phi[j - 1] + d_phi;
  d_d_phi = differentiation_matrix(nphi, 0.0, 2 * pi / nfp);

  // Initialize the axis shape.
  R0 = R0c[0];
  Z0 = Z0c[0];
  R0p = 0.0;
  Z0p = 0.0;
  R0pp = 0.0;
  Z0pp = 0.0;
  R0ppp = 0.0;
  Z0ppp = 0.0;
  for (n = 1; n < R0c.size(); n++) {
    //sinangle = sin(qscfloat(n * nfp) * phi);
    //cosangle = cos(qscfloat(n * nfp) * phi);
    sinangle = sin((n * nfp) * phi);
    cosangle = cos((n * nfp) * phi);
    R0 += R0c[n] * cosangle + R0s[n] * sinangle;
    Z0 += Z0c[n] * cosangle + Z0s[n] * sinangle;
    //R0_extended = R0_extended + R0c[n] * cosangle + R0s[n] * sinangle;
    //Z0_extended = Z0_extended + R0c[n] * cosangle + R0s[n] * sinangle;
    R0p += R0c[n] * (-n*nfp)*sinangle + R0s[n] * (n*nfp)*cosangle;
    Z0p += Z0c[n] * (-n*nfp)*sinangle + Z0s[n] * (n*nfp)*cosangle;
    R0pp += R0c[n] * (-n*nfp*n*nfp)*cosangle
      + R0s[n] * (-n*nfp*n*nfp)*sinangle;
    Z0pp += Z0c[n] * (-n*nfp*n*nfp)*cosangle
      + Z0s[n] * (-n*nfp*n*nfp)*sinangle;
    R0ppp += R0c[n] * (n*nfp*n*nfp*n*nfp)*sinangle
      + R0s[n] * (-n*nfp*n*nfp*n*nfp)*cosangle;
    Z0ppp += Z0c[n] * (n*nfp*n*nfp*n*nfp)*sinangle
      + Z0s[n] * (-n*nfp*n*nfp*n*nfp)*cosangle;
  }
  d_l_d_phi = sqrt(R0 * R0 + R0p * R0p + Z0p * Z0p);
  d2_l_d_phi2 = (R0 * R0p + R0p * R0pp + Z0p * Z0pp) / d_l_d_phi;
  B0_over_abs_G0 = nphi / d_l_d_phi.sum();
  abs_G0_over_B0 = 1 / B0_over_abs_G0;
  d_l_d_varphi = abs_G0_over_B0;
  G0 = sG * abs_G0_over_B0 * B0;
  if (verbose > 0) std::cout << "G0: " << G0 << std::endl;

  // Define some aliases:
#define d_r_d_phi_cylindrical1 R0p
#define d_r_d_phi_cylindrical2 R0
#define d_r_d_phi_cylindrical3 Z0p

#define d2_r_d_phi2_cylindrical1 (R0pp - R0)
#define d2_r_d_phi2_cylindrical2 (2 * R0p)
#define d2_r_d_phi2_cylindrical3 Z0pp

#define d3_r_d_phi3_cylindrical1 (R0ppp - 3 * R0p)
#define d3_r_d_phi3_cylindrical2 (3 * R0pp - R0)
#define d3_r_d_phi3_cylindrical3 Z0ppp

  tangent_cylindrical1 = d_r_d_phi_cylindrical1 / d_l_d_phi;
  tangent_cylindrical2 = d_r_d_phi_cylindrical2 / d_l_d_phi;
  tangent_cylindrical3 = d_r_d_phi_cylindrical3 / d_l_d_phi;

  d_tangent_d_l_cylindrical1 = (-d_r_d_phi_cylindrical1 * d2_l_d_phi2 / d_l_d_phi 
				+ d2_r_d_phi2_cylindrical1) / (d_l_d_phi * d_l_d_phi);
  
  d_tangent_d_l_cylindrical2 = (-d_r_d_phi_cylindrical2 * d2_l_d_phi2 / d_l_d_phi 
				+ d2_r_d_phi2_cylindrical2) / (d_l_d_phi * d_l_d_phi);
  
  d_tangent_d_l_cylindrical3 = (-d_r_d_phi_cylindrical3 * d2_l_d_phi2 / d_l_d_phi 
				+ d2_r_d_phi2_cylindrical3) / (d_l_d_phi * d_l_d_phi);
  /*
  std::cout << "R0: " << R0 << std::endl;
  std::cout << "Z0: " << Z0 << std::endl;
  std::cout << "R0p: " << R0p << std::endl;
  std::cout << "Z0p: " << Z0p << std::endl;
  std::cout << "R0pp: " << R0pp << std::endl;
  std::cout << "Z0pp: " << Z0pp << std::endl;
  std::cout << "R0ppp: " << R0ppp << std::endl;
  std::cout << "Z0ppp: " << Z0ppp << std::endl;
  std::cout << "d_r_d_phi_cylindrical1: " << d_r_d_phi_cylindrical1 << std::endl;
  std::cout << "d_r_d_phi_cylindrical2: " << d_r_d_phi_cylindrical2 << std::endl;
  std::cout << "d_r_d_phi_cylindrical3: " << d_r_d_phi_cylindrical3 << std::endl;
  std::cout << "d_tangent_d_l_cylindrical1: " << d_tangent_d_l_cylindrical1 << std::endl;
  std::cout << "d_tangent_d_l_cylindrical2: " << d_tangent_d_l_cylindrical2 << std::endl;
  std::cout << "d_tangent_d_l_cylindrical3: " << d_tangent_d_l_cylindrical3 << std::endl;
  */
  curvature = sqrt(d_tangent_d_l_cylindrical1 * d_tangent_d_l_cylindrical1 +
		   d_tangent_d_l_cylindrical2 * d_tangent_d_l_cylindrical2 +
		   d_tangent_d_l_cylindrical3 * d_tangent_d_l_cylindrical3);

  axis_length = d_l_d_phi.sum() * d_phi * nfp;
  grid_max_curvature = curvature.max();
  tempvec = curvature * curvature * d_l_d_phi;
  rms_curvature = sqrt((tempvec.sum() * d_phi * nfp) / axis_length);
  // At this point in the fortran code, we find the exact max curvature.
  
  tempvec = R0 * d_l_d_phi;
  mean_R = tempvec.sum() * d_phi * nfp / axis_length;
  tempvec = Z0 * d_l_d_phi;
  mean_Z = tempvec.sum() * d_phi * nfp / axis_length;
  tempvec = (R0 - mean_R) * (R0 - mean_R) * d_l_d_phi;
  standard_deviation_of_R = sqrt(tempvec.sum() * d_phi * nfp / axis_length);
  tempvec = (Z0 - mean_Z) * (Z0 - mean_Z) * d_l_d_phi;
  standard_deviation_of_Z = sqrt(tempvec.sum() * d_phi * nfp / axis_length);

  normal_cylindrical1 = d_tangent_d_l_cylindrical1 / curvature;
  normal_cylindrical2 = d_tangent_d_l_cylindrical2 / curvature;
  normal_cylindrical3 = d_tangent_d_l_cylindrical3 / curvature;
  
  // Compute the vector b = t x n
  binormal_cylindrical1 = tangent_cylindrical2 * normal_cylindrical3 - tangent_cylindrical3 * normal_cylindrical2;
  binormal_cylindrical2 = tangent_cylindrical3 * normal_cylindrical1 - tangent_cylindrical1 * normal_cylindrical3;
  binormal_cylindrical3 = tangent_cylindrical1 * normal_cylindrical2 - tangent_cylindrical2 * normal_cylindrical1;

  // We use the same sign convention for torsion as the Landreman-Sengupta-Plunk paper, wikipedia, and mathworld.wolfram.com/Torsion.html.
  // This sign convention is opposite to Garren & Boozer's sign convention!
  torsion_numerator = (0
		       + d_r_d_phi_cylindrical1 * (d2_r_d_phi2_cylindrical2 * d3_r_d_phi3_cylindrical3 - d2_r_d_phi2_cylindrical3 * d3_r_d_phi3_cylindrical2) 
		       + d_r_d_phi_cylindrical2 * (d2_r_d_phi2_cylindrical3 * d3_r_d_phi3_cylindrical1 - d2_r_d_phi2_cylindrical1 * d3_r_d_phi3_cylindrical3) 
		       + d_r_d_phi_cylindrical3 * (d2_r_d_phi2_cylindrical1 * d3_r_d_phi3_cylindrical2 - d2_r_d_phi2_cylindrical2 * d3_r_d_phi3_cylindrical1));

  tempvec1 = d_r_d_phi_cylindrical2 * d2_r_d_phi2_cylindrical3 - d_r_d_phi_cylindrical3 * d2_r_d_phi2_cylindrical2;
  tempvec2 = d_r_d_phi_cylindrical3 * d2_r_d_phi2_cylindrical1 - d_r_d_phi_cylindrical1 * d2_r_d_phi2_cylindrical3;
  tempvec3 = d_r_d_phi_cylindrical1 * d2_r_d_phi2_cylindrical2 - d_r_d_phi_cylindrical2 * d2_r_d_phi2_cylindrical1;
  torsion_denominator = tempvec1 * tempvec1 + tempvec2 * tempvec2 + tempvec3 * tempvec3;
  
  torsion = torsion_numerator / torsion_denominator;

  etabar_squared_over_curvature_squared = (eta_bar * eta_bar) / (curvature * curvature);
  for (k = 0; k < nphi; k++) {
    for (j = 0; j < nphi; j++) {
      d_d_varphi(j, k) = d_d_phi(j, k) / (B0_over_abs_G0 * d_l_d_phi[j]);
    }
  }

  X1s = 0;
  X1c = eta_bar / curvature;

  // Compute the Boozer toroidal angle along the axis, which is
  // proportional (for QA or QH) to arclength along the axis.
  Boozer_toroidal_angle[0] = 0.0;
  for (j = 1; j < nphi; j++) {
    // To get toroidal angle on the full mesh, we need d_l_d_phi on
    // the half mesh.
    Boozer_toroidal_angle[j] = Boozer_toroidal_angle[j - 1]
      + (d_l_d_phi[j - 1] + d_l_d_phi[j]);
  }
  Boozer_toroidal_angle *= (0.5 * d_phi * 2 * pi / axis_length);

  calculate_helicity();

  if (verbose > 0) {
    std::cout << "R0c:" << R0c << std::endl;
    std::cout << "R0s:" << R0s << std::endl;
    std::cout << "Z0c:" << Z0c << std::endl;
    std::cout << "Z0s:" << Z0s << std::endl;
    std::cout << "phi:" << phi << std::endl << std::endl;
    std::cout << "R0:" << R0 << std::endl << std::endl;
    std::cout << "curvature:" << curvature << std::endl << std::endl;
    std::cout << "torsion:" << torsion << std::endl << std::endl;
    std::cout << "stddevR: " << standard_deviation_of_R
	      << "  stddevZ: " << standard_deviation_of_Z << std::endl;
  }

  if (verbose > 0) {
    end_time = std::clock();
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for init_axis from chrono:           "
              << elapsed.count() << " seconds" << std::endl;
    std::cout << "Time for init_axis from ctime (CPU time): "
              << double(end_time - start_time) / CLOCKS_PER_SEC
              << " seconds" << std::endl;
  }

}
