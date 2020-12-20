#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

/** This test case compares the grad B tensor to eq (3.14) in
    Landreman JPP (2021) for an axisymmetric vacuum field.
 */
TEST_CASE("grad B tensor for an axisymmetric vacuum field") {
  Qsc q;
  qscfloat val;
  q.I2 = 0.0;
  q.sigma0 = 0.0;
  q.verbose = 0;

  q.R0c.resize(1, 1.0);
  q.R0s.resize(1, 0.0);
  q.Z0c.resize(1, 0.0);
  q.Z0s.resize(1, 0.0);

  /*
  for (int sG = -1; sG <= -1; sG += 2) {
    for (int spsi = -1; spsi <= -1; spsi += 2) {
      for (int nfp = 1; nfp < 2; nfp++) {
	for (qscfloat eta_bar = 0.6; eta_bar < 0.61; eta_bar += 0.31) {
	  for (int nphi = 3; nphi < 4; nphi += 2) {
	    for (qscfloat B0 = 0.7; B0 < 0.71; B0 += 0.32) {
	      for (qscfloat R0 = 0.65; R0 < 0.66; R0 += 0.34) {
  */
  
  for (int sG = -1; sG <= 1; sG += 2) {
    for (int spsi = -1; spsi <= 1; spsi += 2) {
      for (int nfp = 1; nfp < 4; nfp++) {
	for (qscfloat eta_bar = 0.6; eta_bar < 2.0; eta_bar += 0.31) {
	  for (int nphi = 3; nphi < 10; nphi += 2) {
	    for (qscfloat B0 = 0.7; B0 < 1.8; B0 += 0.32) {
	      for (qscfloat R0 = 0.65; R0 < 1.6; R0 += 0.34) {
		q.sG = sG;
		q.spsi = spsi;
		q.nfp = nfp;
		q.eta_bar = eta_bar;
		q.nphi = nphi;
		q.R0c[0] = R0;
		q.B0 = B0;
		
		q.calculate();

		val = sG * B0 / R0;
		for (int j = 0; j < nphi; j++) {
		  // In the dimensions of size 3, the order of elements is
		  // (normal, binormal, tangent)
		  
		  // tn
		  CHECK(Approx(q.grad_B_tensor(j, 2, 0)) == val);
		  // nt
		  CHECK(Approx(q.grad_B_tensor(j, 0, 2)) == val);
		  // bb
		  CHECK(Approx(q.grad_B_tensor(j, 1, 1)) == 0.0);
		  // nn
		  CHECK(Approx(q.grad_B_tensor(j, 0, 0)) == 0.0);
		  // bn
		  CHECK(Approx(q.grad_B_tensor(j, 1, 0)) == 0.0);
		  // nb
		  CHECK(Approx(q.grad_B_tensor(j, 0, 1)) == 0.0);
		  // tt
		  CHECK(Approx(q.grad_B_tensor(j, 2, 2)) == 0.0);
		  // tb
		  CHECK(Approx(q.grad_B_tensor(j, 2, 1)) == 0.0);
		  // bt
		  CHECK(Approx(q.grad_B_tensor(j, 1, 2)) == 0.0);

		  CHECK(Approx(q.L_grad_B[j]) == R0);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

