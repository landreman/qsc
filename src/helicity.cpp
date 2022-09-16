#include <stdexcept>
#include "qsc.hpp"

using namespace qsc;

void Qsc::calculate_helicity() {
  int j;
  qscfloat normal_cylindrical_R;
  for (j = 0; j < nphi; j++) {
    normal_cylindrical_R = normal_cartesian1[j] * cosphi[j] + normal_cartesian2[j] * sinphi[j];
    if (normal_cylindrical_R >= 0) {
      if (normal_cartesian3[j] >= 0) {
	quadrant[j] = 1;
      } else {
	quadrant[j] = 4;
      }
    } else {
      if (normal_cartesian3[j] >= 0) {
	quadrant[j] = 2;
      } else {
	quadrant[j] = 3;
      }
    }
  }
  quadrant[nphi] = quadrant[0];

  int counter = 0;
  for (j = 0; j < nphi; j++) {
    if (quadrant[j] == 4 && quadrant[j + 1] == 1) {
      counter++;
    } else if (quadrant[j] == 1 && quadrant[j + 1] == 4) {
      counter--;
    } else {
      counter += quadrant[j + 1] - quadrant[j];
    }
  }

  // It is necessary to flip the sign of axis_helicity in order to
  // maintain "iota_N = iota + axis_helicity" under the parity
  // transformations.
  counter *= spsi * sG;

  helicity = counter / 4;
  if (verbose > 0) {
    std::cout << "Axis helicity counter: " << counter << "  helicity: " << helicity << std::endl;
  }
  if (counter % 4 != 0) {
    throw std::runtime_error("Axis helicity counter was not a multiple of 4");
  }
}
