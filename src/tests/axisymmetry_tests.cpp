#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

/** See the document 20190214-01 Setting up Shafranov shift
 *  benchmark.docx and 20190214-01 Setting up Shafranov shift
 *  benchmark.nb to see how to set the constants below.  We can pick
 *  R0c, B0, and I2 freely. p2 is free, but we want it close in
 *  magnitude to the term in the Shafranov shift that is independent
 *  of pressure. Then eta_bar and B2c are fixed.
*/
TEST_CASE("Compare Shafranov shift to VMEC and textbook expressions for a circular-cross-section tokamak") {
  Qsc q;

  q.nfp = 1;
  q.nphi = 3;
  q.order_r_option = "r2";

  q.R0c.resize(1, 0.0);
  q.R0s.resize(1, 0.0);
  q.Z0c.resize(1, 0.0);
  q.Z0s.resize(1, 0.0);

  qscfloat R0 = 1.7;
  q.R0c[0] = R0;
  q.B0 = 0.8;
  q.I2 = 1.04;
  // For the p2 and p-independent terms in Shafranov shift to be comparable, we should choose
  // p2 ~ -I^2 / (4 mu0)
  q.p2 = -3.0e+5;
  
  q.eta_bar = 1 / R0; // So elongation is 1.
  q.B2c = mu0 * q.p2 * q.B0 / (2 * q.I2 * q.I2 * R0 * R0) - q.B0 / (8 * R0 * R0);

  for (q.sG = -1; q.sG <= 1; q.sG += 2) {
    for (q.spsi = -1; q.spsi <= 1; q.spsi += 2) {
      CAPTURE(q.sG);
      CAPTURE(q.spsi);
      q.init();
      q.calculate();

      CHECK(Approx(q.iota) == q.sG * q.I2 * R0 / q.B0);
      for (int j = 0; j < q.nphi; j++) {
	CAPTURE(j);
	CHECK(Approx(q.X2c[j]) == mu0 * q.p2 / (2 * q.I2 * q.I2 * R0) - 5 / (8 * R0));
	CHECK(Approx(q.X20[j]) == -mu0 * q.p2 / (q.I2 * q.I2 * R0) + 3 / (4 * R0));
	CHECK(Approx(q.B20[j]) == q.B0 * (-mu0 * q.p2 / (q.B0 * q.B0)
					  - mu0 * q.p2 / (q.I2 * q.I2 * R0 * R0)
					  + 5 / (4 * R0 * R0)
					  - q.I2 * q.I2 / (2 * q.B0 * q.B0)));
	CHECK(Approx(q.Y2s[j]) == q.sG * q.spsi * (mu0 * q.p2 / (2 * q.I2 * q.I2 * R0) - 5 / (8 * R0)));
      }
    }
  }
}

