#include <chrono>
#include "qsc.hpp"

using namespace qsc;

/** Allocate all of the arrays (Vectors), matrices, and higher-rank tensors that will be used.
 */
void Qsc::allocate() {
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) start = std::chrono::steady_clock::now();
  
  // Ensure nphi is odd:
  if (nphi % 2 == 0) nphi++;

  phi0.resize(nphi, 0.0);
  phi.resize(nphi, 0.0);
  sinphi.resize(nphi, 0.0);
  cosphi.resize(nphi, 0.0);
  R0.resize(nphi, 0.0);
  Z0.resize(nphi, 0.0);
  R0p.resize(nphi, 0.0);
  Z0p.resize(nphi, 0.0);
  R0pp.resize(nphi, 0.0);
  Z0pp.resize(nphi, 0.0);
  R0ppp.resize(nphi, 0.0);
  Z0ppp.resize(nphi, 0.0);
  f.resize(nphi, 0.0);
  fp.resize(nphi, 0.0);
  fpp.resize(nphi, 0.0);
  fppp.resize(nphi, 0.0);
  curvature.resize(nphi, 0.0);
  torsion.resize(nphi, 0.0);
  sinangle.resize(nphi, 0.0);
  cosangle.resize(nphi, 0.0);
  d_l_d_phi0.resize(nphi, 0.0);
  d2_l_d_phi02.resize(nphi, 0.0);
  
  tangent_cartesian1.resize(nphi, 0.0);
  tangent_cartesian2.resize(nphi, 0.0);
  tangent_cartesian3.resize(nphi, 0.0);

  normal_cartesian1.resize(nphi, 0.0);
  normal_cartesian2.resize(nphi, 0.0);
  normal_cartesian3.resize(nphi, 0.0);

  binormal_cartesian1.resize(nphi, 0.0);
  binormal_cartesian2.resize(nphi, 0.0);
  binormal_cartesian3.resize(nphi, 0.0);

  d_tangent_d_l_cartesian1.resize(nphi, 0.0);
  d_tangent_d_l_cartesian2.resize(nphi, 0.0);
  d_tangent_d_l_cartesian3.resize(nphi, 0.0);

  tempvec.resize(nphi, 0.0);
  tempvec1.resize(nphi, 0.0);
  tempvec2.resize(nphi, 0.0);
  tempvec3.resize(nphi, 0.0);

  torsion_numerator.resize(nphi, 0.0);
  torsion_denominator.resize(nphi, 0.0);

  Boozer_toroidal_angle.resize(nphi, 0.0);
  etabar_squared_over_curvature_squared.resize(nphi, 0.0);

  d_d_phi0.resize(nphi, nphi, 0.0);
  d_d_varphi.resize(nphi, nphi, 0.0);
  work_matrix.resize(nphi, nphi, 0.0);
  
  X1s.resize(nphi, 0.0);
  X1c.resize(nphi, 0.0);
  Y1s.resize(nphi, 0.0);
  Y1c.resize(nphi, 0.0);
  sigma.resize(nphi, 0.0);
  elongation.resize(nphi, 0.0);
  
  quadrant.resize(nphi + 1, 0);
  state.resize(nphi, 0);
  residual.resize(nphi, 0);
  work1.resize(nphi, 0.0);
  work2.resize(nphi, 0.0);
  ipiv.resize(nphi, 0);

  d_X1c_d_varphi.resize(nphi, 0.0);
  d_Y1s_d_varphi.resize(nphi, 0.0);
  d_Y1c_d_varphi.resize(nphi, 0.0);

  grad_B_tensor.resize(nphi, 3, 3, 0.0);
  L_grad_B.resize(nphi, 0.0);
  L_grad_B_inverse.resize(nphi, 0.0);

  if (at_least_order_r2) {
    B20.resize(nphi, 0.0);
    B20_anomaly.resize(nphi, 0.0);

    X20.resize(nphi, 0.0);
    X2s.resize(nphi, 0.0);
    X2c.resize(nphi, 0.0);

    Y20.resize(nphi, 0.0);
    Y2s.resize(nphi, 0.0);
    Y2c.resize(nphi, 0.0);

    Z20.resize(nphi, 0.0);
    Z2s.resize(nphi, 0.0);
    Z2c.resize(nphi, 0.0);

    V1.resize(nphi, 0.0);
    V2.resize(nphi, 0.0);
    V3.resize(nphi, 0.0);
    
    rs.resize(nphi, 0.0);
    rc.resize(nphi, 0.0);
    qs.resize(nphi, 0.0);
    qc.resize(nphi, 0.0);

    r2_matrix.resize(2 * nphi, 2 * nphi, 0.0);
    r2_rhs.resize(2 * nphi, 0.0);
    r2_ipiv.resize(2 * nphi, 0);

    Y2s_from_X20.resize(nphi, 0.0);
    Y2s_inhomogeneous.resize(nphi, 0.0);
    Y2c_from_X20.resize(nphi, 0.0);
    Y2c_inhomogeneous.resize(nphi, 0.0);

    fX0_from_X20.resize(nphi, 0.0);
    fX0_from_Y20.resize(nphi, 0.0);
    fX0_inhomogeneous.resize(nphi, 0.0);

    fXs_from_X20.resize(nphi, 0.0);
    fXs_from_Y20.resize(nphi, 0.0);
    fXs_inhomogeneous.resize(nphi, 0.0);

    fXc_from_X20.resize(nphi, 0.0);
    fXc_from_Y20.resize(nphi, 0.0);
    fXc_inhomogeneous.resize(nphi, 0.0);

    fY0_from_X20.resize(nphi, 0.0);
    fY0_from_Y20.resize(nphi, 0.0);
    fY0_inhomogeneous.resize(nphi, 0.0);

    fYs_from_X20.resize(nphi, 0.0);
    fYs_from_Y20.resize(nphi, 0.0);
    fYs_inhomogeneous.resize(nphi, 0.0);

    fYc_from_X20.resize(nphi, 0.0);
    fYc_from_Y20.resize(nphi, 0.0);
    fYc_inhomogeneous.resize(nphi, 0.0);

    d_curvature_d_varphi.resize(nphi, 0.0);
    d_torsion_d_varphi.resize(nphi, 0.0);
    
    d_X20_d_varphi.resize(nphi, 0.0);
    d_X2s_d_varphi.resize(nphi, 0.0);
    d_X2c_d_varphi.resize(nphi, 0.0);
    
    d_Y20_d_varphi.resize(nphi, 0.0);
    d_Y2s_d_varphi.resize(nphi, 0.0);
    d_Y2c_d_varphi.resize(nphi, 0.0);
    
    d_Z20_d_varphi.resize(nphi, 0.0);
    d_Z2s_d_varphi.resize(nphi, 0.0);
    d_Z2c_d_varphi.resize(nphi, 0.0);
    
    d2_X1c_d_varphi2.resize(nphi, 0.0);
    d2_Y1c_d_varphi2.resize(nphi, 0.0);
    d2_Y1s_d_varphi2.resize(nphi, 0.0);
    
    d2_X20_d_varphi2.resize(nphi, 0.0);
    d2_X2s_d_varphi2.resize(nphi, 0.0);
    d2_X2c_d_varphi2.resize(nphi, 0.0);
    
    d2_Y20_d_varphi2.resize(nphi, 0.0);
    d2_Y2s_d_varphi2.resize(nphi, 0.0);
    d2_Y2c_d_varphi2.resize(nphi, 0.0);
    
    d2_Z20_d_varphi2.resize(nphi, 0.0);
    d2_Z2s_d_varphi2.resize(nphi, 0.0);
    d2_Z2c_d_varphi2.resize(nphi, 0.0);

    grad_grad_B_tensor.resize(nphi, 3, 3, 3, 0.0);
    L_grad_grad_B.resize(nphi, 0.0);
    L_grad_grad_B_inverse.resize(nphi, 0.0);

    r_hat_singularity_robust.resize(nphi, 0.0);
  }

  if (order_r2p1) {
    lambda_for_XY3.resize(nphi, 0.0);
  }
  
  if (order_r3) {
    X3c1.resize(nphi, 0.0);
    X3c3.resize(nphi, 0.0);
    X3s1.resize(nphi, 0.0);
    X3s3.resize(nphi, 0.0);
    
    Y3c1.resize(nphi, 0.0);
    Y3c3.resize(nphi, 0.0);
    Y3s1.resize(nphi, 0.0);
    Y3s3.resize(nphi, 0.0);

    Z3c1.resize(nphi, 0.0);
    Z3c3.resize(nphi, 0.0);
    Z3s1.resize(nphi, 0.0);
    Z3s3.resize(nphi, 0.0);

    d_X3c1_d_varphi.resize(nphi, 0.0);
    d_X3c3_d_varphi.resize(nphi, 0.0);
    d_X3s1_d_varphi.resize(nphi, 0.0);
    d_X3s3_d_varphi.resize(nphi, 0.0);

    d_Y3c1_d_varphi.resize(nphi, 0.0);
    d_Y3c3_d_varphi.resize(nphi, 0.0);
    d_Y3s1_d_varphi.resize(nphi, 0.0);
    d_Y3s3_d_varphi.resize(nphi, 0.0);

    d2_X3c1_d_varphi2.resize(nphi, 0.0);
    d2_X3c3_d_varphi2.resize(nphi, 0.0);
    d2_X3s1_d_varphi2.resize(nphi, 0.0);
    d2_X3s3_d_varphi2.resize(nphi, 0.0);

    d2_Y3c1_d_varphi2.resize(nphi, 0.0);
    d2_Y3c3_d_varphi2.resize(nphi, 0.0);
    d2_Y3s1_d_varphi2.resize(nphi, 0.0);
    d2_Y3s3_d_varphi2.resize(nphi, 0.0);
}
  
  if (verbose > 0) {
    auto end = std::chrono::steady_clock::now();    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for allocate: "
              << elapsed.count() << " seconds" << std::endl;
  }
}
