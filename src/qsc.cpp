#include "qsc.hpp"

using namespace qsc;

Qsc::Qsc() :
  // Call constructor of member objects:
  d_d_phi(1, 1),
  d_d_varphi(1, 1)
{
  // Set defaults.
  verbose = 1;
  
  sG = 1;
  spsi = 1;
  B0 = 1.0;
  eta_bar = -1.0;

  nfp = 3;
  nphi = 15;
  
  R0c.resize(1, 1.0);
  R0s.resize(1, 0.0);
  Z0c.resize(1, 0.0);
  Z0s.resize(1, 0.0);

}
