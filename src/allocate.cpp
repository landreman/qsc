#include "qsc.hpp"

using namespace qsc;

/** Allocate all of the arrays (Vectors) that will be used.
 */
void Qsc::allocate() {
  // Ensure nphi is odd:
  if (nphi % 2 == 0) nphi++;

  phi.resize(nphi, 0.0);
  R0.resize(nphi, 0.0);
  Z0.resize(nphi, 0.0);
  R0p.resize(nphi, 0.0);
  Z0p.resize(nphi, 0.0);
  R0pp.resize(nphi, 0.0);
  Z0pp.resize(nphi, 0.0);
  R0ppp.resize(nphi, 0.0);
  Z0ppp.resize(nphi, 0.0);
  curvature.resize(nphi, 0.0);
  torsion.resize(nphi, 0.0);
  sinangle.resize(nphi, 0.0);
  cosangle.resize(nphi, 0.0);
  d_l_d_phi.resize(nphi, 0.0);
  d2_l_d_phi2.resize(nphi, 0.0);
  
  tangent_cylindrical1.resize(nphi, 0.0);
  tangent_cylindrical2.resize(nphi, 0.0);
  tangent_cylindrical3.resize(nphi, 0.0);

  normal_cylindrical1.resize(nphi, 0.0);
  normal_cylindrical2.resize(nphi, 0.0);
  normal_cylindrical3.resize(nphi, 0.0);

  binormal_cylindrical1.resize(nphi, 0.0);
  binormal_cylindrical2.resize(nphi, 0.0);
  binormal_cylindrical3.resize(nphi, 0.0);

  d_tangent_d_l_cylindrical1.resize(nphi, 0.0);
  d_tangent_d_l_cylindrical2.resize(nphi, 0.0);
  d_tangent_d_l_cylindrical3.resize(nphi, 0.0);

  tempvec.resize(nphi, 0.0);
  tempvec1.resize(nphi, 0.0);
  tempvec2.resize(nphi, 0.0);
  tempvec3.resize(nphi, 0.0);

  torsion_numerator.resize(nphi, 0.0);
  torsion_denominator.resize(nphi, 0.0);

  Boozer_toroidal_angle.resize(nphi, 0.0);
  B1Squared_over_curvatureSquared.resize(nphi, 0.0);

  d_d_varphi.resize(nphi, nphi, 0.0);
  
  X1s.resize(nphi, 0.0);
  X1c.resize(nphi, 0.0);
  quadrant.resize(nphi + 1, 0);
}