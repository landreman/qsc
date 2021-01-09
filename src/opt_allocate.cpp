#include "opt.hpp"

using namespace qsc;

/** Allocate all of the arrays (Vectors) and matrices that will be used.
 */
void Opt::allocate() {
  q.allocate();
  sqrt_d_l_d_phi.resize(q.nphi, 0.0);

  iter_objective_function.resize(max_iter, 0.0);
  iter_B20_term.resize(max_iter, 0.0);
  iter_iota_term.resize(max_iter, 0.0);
  iter_R0_term.resize(max_iter, 0.0);
  iter_d2_volume_d_psi2_term.resize(max_iter, 0.0);
  iter_XY2_term.resize(max_iter, 0.0);
  iter_XY2Prime_term.resize(max_iter, 0.0);
  iter_XY3_term.resize(max_iter, 0.0);
  iter_XY3Prime_term.resize(max_iter, 0.0);
  iter_grad_grad_B_term.resize(max_iter, 0.0);
  
  iter_eta_bar.resize(max_iter, 0.0);
  iter_sigma0.resize(max_iter, 0.0);
  iter_B2s.resize(max_iter, 0.0);
  iter_B2c.resize(max_iter, 0.0);
  
  iter_R0c.resize(q.R0c.size(), max_iter, 0.0);
  iter_R0s.resize(q.R0c.size(), max_iter, 0.0);
  iter_Z0c.resize(q.R0c.size(), max_iter, 0.0);
  iter_Z0s.resize(q.R0c.size(), max_iter, 0.0);
  
  iter_min_R0.resize(max_iter, 0.0);
  iter_max_curvature.resize(max_iter, 0.0);
  iter_iota.resize(max_iter, 0.0);
  iter_max_elongation.resize(max_iter, 0.0);
  iter_min_L_grad_B.resize(max_iter, 0.0);
  iter_min_L_grad_grad_B.resize(max_iter, 0.0);
  iter_r_singularity.resize(max_iter, 0.0);
  iter_B20_variation.resize(max_iter, 0.0);
  iter_B20_residual.resize(max_iter, 0.0);
  iter_d2_volume_d_psi2.resize(max_iter, 0.0);
  iter_DMerc_times_r2.resize(max_iter, 0.0);
  iter_standard_deviation_of_R.resize(max_iter, 0.0);
  iter_standard_deviation_of_Z.resize(max_iter, 0.0);
  
}
