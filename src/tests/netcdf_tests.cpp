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
  CHECK(Approx(qfile.eta_bar) == q11.eta_bar);
  CHECK(Approx(qfile.B2c) == q11.B2c);
  CHECK(Approx(qfile.B2s) == q11.B2s);
  CHECK(Approx(qfile.p2) == q11.p2);
  CHECK(Approx(qfile.d_phi) == q11.d_phi);
  CHECK(Approx(qfile.B0) == q11.B0);
  CHECK(Approx(qfile.G0) == q11.G0);
  CHECK(qfile.sG == q11.sG);
  CHECK(qfile.spsi == q11.spsi);
  CHECK(Approx(qfile.axis_length) == q11.axis_length);
  CHECK(Approx(qfile.d_l_d_varphi) == q11.d_l_d_varphi);
  CHECK(Approx(qfile.B0_over_abs_G0) == q11.B0_over_abs_G0);
  CHECK(Approx(qfile.abs_G0_over_B0) == q11.abs_G0_over_B0);
  CHECK(Approx(qfile.rms_curvature) == q11.rms_curvature);
  CHECK(Approx(qfile.mean_elongation) == q11.mean_elongation);
  CHECK(Approx(qfile.mean_R) == q11.mean_R);
  CHECK(Approx(qfile.mean_Z) == q11.mean_Z);
  CHECK(Approx(qfile.standard_deviation_of_R) == q11.standard_deviation_of_R);
  CHECK(Approx(qfile.standard_deviation_of_Z) == q11.standard_deviation_of_Z);
  CHECK(Approx(qfile.iota) == q11.iota);
  CHECK(Approx(qfile.iota_N) == q11.iota_N);
  CHECK(Approx(qfile.I2) == q11.I2);
  CHECK(Approx(qfile.sigma0) == q11.sigma0);
  CHECK(Approx(qfile.grid_max_curvature) == q11.grid_max_curvature);
  CHECK(Approx(qfile.grid_max_elongation) == q11.grid_max_elongation);
  CHECK(Approx(qfile.grid_min_R0) == q11.grid_min_R0);
  CHECK(Approx(qfile.newton_tolerance) == q11.newton_tolerance);
  CHECK(qfile.max_newton_iterations == q11.max_newton_iterations);
  CHECK(qfile.max_linesearch_iterations == q11.max_linesearch_iterations);
  CHECK(qfile.helicity == q11.helicity);
  // CHECK(Approx(qfile.) == q11.);

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

