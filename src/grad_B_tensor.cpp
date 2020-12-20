#include "qsc.hpp"

using namespace qsc;

/** Compute the grad B tensor
 */
void Qsc::calculate_grad_B_tensor() {
  // In the dimensions of size 3, the order of elements is
  // (normal, binormal, tangent)
  
  qscfloat factor = spsi * B0 / d_l_d_varphi;

  // Evaluate eq (3.12) in Landreman JPP (2021):
  
  // tn
  tempvec = (sG * B0) * curvature;
  if (verbose > 0) std::cout << "in grad B tensor: tn=" << tempvec << std::endl;
  grad_B_tensor.set_row(tempvec, 2, 0);

  // nt - same as tn
  grad_B_tensor.set_row(tempvec, 0, 2);

  // bb
  tempvec = factor * (X1c * d_Y1s_d_varphi - iota_N * X1c * Y1c);
  grad_B_tensor.set_row(tempvec, 1, 1);

  // nn
  tempvec = factor * (d_X1c_d_varphi * Y1s + iota_N * X1c * Y1c);
  grad_B_tensor.set_row(tempvec, 0, 0);

  // bn
  tempvec = factor * ((-sG * spsi * d_l_d_varphi) * torsion
		      - iota_N * X1c * X1c);
  grad_B_tensor.set_row(tempvec, 1, 0);

  // nb
  tempvec = factor * (d_Y1c_d_varphi * Y1s - d_Y1s_d_varphi * Y1c
		      + (sG * spsi * d_l_d_varphi) * torsion
		      + iota_N * (Y1s * Y1s + Y1c * Y1c));
  grad_B_tensor.set_row(tempvec, 0, 1);

  // Evaluate eq (3.1) in Landreman JPP (2021):

  for (int j = 0; j < nphi; j++) {
    L_grad_B[j] = B0 * sqrt(2 / (grad_B_tensor(j, 2, 0) * grad_B_tensor(j, 2, 0) +
				 grad_B_tensor(j, 0, 2) * grad_B_tensor(j, 0, 2) +
				 grad_B_tensor(j, 1, 1) * grad_B_tensor(j, 1, 1) +
				 grad_B_tensor(j, 0, 0) * grad_B_tensor(j, 0, 0) +
				 grad_B_tensor(j, 1, 0) * grad_B_tensor(j, 1, 0) +
				 grad_B_tensor(j, 0, 1) * grad_B_tensor(j, 0, 1)));
  }
  grid_min_L_grad_B = L_grad_B.min();
}
