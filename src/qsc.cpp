#include "qsc.hpp"

using namespace qsc;

Qsc::Qsc() {
  // Set defaults.
  sG = 1;
  spsi = 1;
  B0 = 1.0;

  R0c.resize(1, 1.0);
  R0s.resize(1, 0.0);
  Z0c.resize(1, 0.0);
  Z0s.resize(1, 0.0);
}
