#include <string>
#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

/** Save Qsc data to a NetCDF file. Read the data in to a different
    Qsc object. Verify that the two objects now have the same data.
 */
TEST_CASE("netcdf") {
  Qsc q11("r1 section 5.1");
  Qsc qfile;
  std::string filename = "qsc_out.unitTests.nc";
  q11.write_netcdf(filename);
  qfile.read_netcdf(filename, 'C');

  // Scalars
  CHECK(qfile.nphi == q11.nphi);
  CHECK(qfile.nfp == q11.nfp);
  CHECK(qfile.helicity == q11.helicity);
  CHECK(Approx(qfile.iota) == q11.iota);
  CHECK(Approx(qfile.mean_elongation) == q11.mean_elongation);
  CHECK(Approx(qfile.rms_curvature) == q11.rms_curvature);

  // Vectors
  for (int j = 0; j < q11.nphi; j++) {
    CHECK(Approx(qfile.phi[j]) == q11.phi[j]);
    CHECK(Approx(qfile.curvature[j]) == q11.curvature[j]);
    CHECK(Approx(qfile.torsion[j]) == q11.torsion[j]);
    CHECK(Approx(qfile.sigma[j]) == q11.sigma[j]);
    CHECK(Approx(qfile.X1c[j]) == q11.X1c[j]);
    CHECK(Approx(qfile.Y1c[j]) == q11.Y1c[j]);
    CHECK(Approx(qfile.Y1s[j]) == q11.Y1s[j]);
    CHECK(Approx(qfile.R0[j]) == q11.R0[j]);
    CHECK(Approx(qfile.Z0[j]) == q11.Z0[j]);
    CHECK(Approx(qfile.d_l_d_phi[j]) == q11.d_l_d_phi[j]);
    CHECK(Approx(qfile.d2_l_d_phi2[j]) == q11.d2_l_d_phi2[j]);
    CHECK(Approx(qfile.elongation[j]) == q11.elongation[j]);
    CHECK(Approx(qfile.Boozer_toroidal_angle[j]) == q11.Boozer_toroidal_angle[j]);
    CHECK(Approx(qfile.d_X1c_d_varphi[j]) == q11.d_X1c_d_varphi[j]);
    CHECK(Approx(qfile.d_Y1c_d_varphi[j]) == q11.d_Y1c_d_varphi[j]);
    CHECK(Approx(qfile.d_Y1s_d_varphi[j]) == q11.d_Y1s_d_varphi[j]);
    CHECK(Approx(qfile.L_grad_B[j]) == q11.L_grad_B[j]);
    CHECK(Approx(qfile.L_grad_B_inverse[j]) == q11.L_grad_B_inverse[j]);
  }
}

