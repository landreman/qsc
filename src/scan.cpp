#include "scan.hpp"

using namespace qsc;

void Scan::defaults() {
  // Set defaults.
  max_seconds = 60;

  max_keep_per_proc = 1000000000; // 1e9

  keep_all = true;
  min_iota_to_keep = -1.0;
  min_R0_to_keep = -1.0;
  max_elongation_to_keep = 10.0;
  min_L_grad_B_to_keep = -1.0;
  min_L_grad_grad_B_to_keep = -1.0;
  max_d2_volume_d_psi2_to_keep = 1.0e+30;
  min_DMerc_to_keep = -1.0e+30;
}

Scan::Scan() {
  defaults();
  /*  
  R0c.resize(1, 1.0);
  R0s.resize(1, 0.0);
  Z0c.resize(1, 0.0);
  Z0s.resize(1, 0.0);
  */
}
