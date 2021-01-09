#include <vector>
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
  
  qscfloat tol = 1.0e-12;
  if (single) tol = 2.0e-5;

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
	      for (qscfloat sigma0 = 0; sigma0 < 1; sigma0 += 0.71) {
		for (qscfloat R0 = 0.65; R0 < 1.6; R0 += 0.34) {
		  q.sG = sG;
		  q.spsi = spsi;
		  q.nfp = nfp;
		  q.eta_bar = eta_bar;
		  q.nphi = nphi;
		  q.R0c[0] = R0;
		  q.B0 = B0;
		  q.sigma0 = sigma0;
		  q.max_linesearch_iterations = 20; // It helps to use this increased value in single precision

		  q.init();
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
		    CHECK(Approx(q.grad_B_tensor(j, 0, 1)).epsilon(tol) == 0.0); // For some reason, the error in this term in single precision is a bit larger, so we need a wider tolerance than the default.
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
}

/** This test case compares the grad grad B tensor to eq (3.15) in
    Landreman JPP (2021) for an axisymmetric vacuum field.
 */
TEST_CASE("grad grad B tensor for an axisymmetric vacuum field") {
  Qsc q;
  qscfloat val, tol, zero_tol;
  // A slight nonzero iota is needed, or else the O(r^2) linear solve is rank-deficient
  if (single) {
    q.I2 = 1.0e-4;
    zero_tol = 1.0e-2;
    tol = 1.0e-4;
  } else {
    q.I2 = 1.0e-7;
    zero_tol = 1.0e-5;
    tol = 1.0e-8;
  }
  q.sigma0 = 0.0;
  q.verbose = 0;
  q.order_r_option = ORDER_R_OPTION_R2;
  
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
      for (int nfp = 1; nfp < 4; nfp += 2) {
	for (qscfloat eta_bar = 0.6; eta_bar < 2.0; eta_bar += 0.31) {
	  for (int nphi = 3; nphi < 10; nphi += 4) {
	    for (qscfloat sigma0 = 0; sigma0 < 1; sigma0 += 0.71) {
	      for (qscfloat B0 = 0.7; B0 < 1.8; B0 += 0.32) {
		for (qscfloat B2c = -0.41; B2c < 1.2; B2c += 1.5) {
		  for (qscfloat B2s = 0.0; B2s < 1.0; B2s += 0.63) {
		    for (qscfloat R0 = 0.65; R0 < 1.6; R0 += 0.34) {
		      q.sG = sG;
		      q.spsi = spsi;
		      q.nfp = nfp;
		      q.eta_bar = eta_bar;
		      q.nphi = nphi;
		      q.R0c[0] = R0;
		      q.B0 = B0;
		      q.sigma0 = sigma0;
		      q.B2c = B2c;
		      q.B2s = B2s;

		      q.init();
		      q.calculate();
		      
		      // Prefactor in eq (3.15):
		      val = 2 * sG * B0 / (R0 * R0);
		      
		      for (int jphi = 0; jphi < nphi; jphi++) {
			CAPTURE(jphi);
			for (int j1 = 0; j1 < 3; j1++) {
			  CAPTURE(j1);
			  for (int j2 = 0; j2 < 3; j2++) {
			    CAPTURE(j2);
			    for (int j3 = 0; j3 < 3; j3++) {
			      CAPTURE(j3);
			      if (j1 == 2 && j2 == 2 && j3 == 2) {
				// ttt component
				CHECK(Approx(q.grad_grad_B_tensor(jphi, j1, j2, j3)).epsilon(tol) == -val);
			      } else if ((j1 == 2 && j2 == 0 && j3 == 0) ||
					 (j1 == 0 && j2 == 2 && j3 == 0) ||
					 (j1 == 0 && j2 == 0 && j3 == 2)) {
				// tnn, ntn, and nnt components
				CHECK(Approx(q.grad_grad_B_tensor(jphi, j1, j2, j3)).epsilon(tol) == val);
			      } else {
				// All other components
				CHECK(Approx(q.grad_grad_B_tensor(jphi, j1, j2, j3)).epsilon(zero_tol) == 0.0);
			      }
			    }
			  }
			}
			CHECK(Approx(q.L_grad_grad_B[jphi]) == R0);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

/** Compare Rogerio's derivation to mine to ensure the results
 * coincide.  Also verify symmetry in the first two components.
 */
TEST_CASE("grad grad B tensor alternative derivation and symmetry") {
  std::vector<std::string> configs = {
    "r2 section 5.1",
    "r2 section 5.2",
    "r2 section 5.3",
    "r2 section 5.4",
    "r2 section 5.5"};

  qscfloat tol;
  if (single) {
    tol = 3.0e-3;
  } else {
    tol = 1.0e-8;
  }
  for (int jconfig = 0; jconfig < configs.size(); jconfig++) {
    std::cout << "jconfig=" << jconfig << std::endl;
    CAPTURE(jconfig);
    
    Qsc q(configs[jconfig]);
    q.nphi = 65;
    q.max_linesearch_iterations = 20; // It helps to have this larger value in single precision
  
    for (int sG = -1; sG <= 1; sG += 2) {
      for (int spsi = -1; spsi <= 1; spsi += 2) {
	for (qscfloat B0 = 0.9; B0 < 2.0; B0 += 0.7) {
	  q.sG = sG;
	  q.spsi = spsi;
	  q.B0 = B0;
	  q.init();
	  q.calculate();
	  auto tensor = q.calculate_grad_grad_B_tensor_alt();
	  
	  for (int j3 = 0; j3 < 3; j3++) {
	    CAPTURE(j3);
	    for (int j2 = 0; j2 < 3; j2++) {
	      CAPTURE(j2);
	      for (int j1 = 0; j1 < 3; j1++) {
		CAPTURE(j1);
		for (int jphi = 0; jphi < q.nphi; jphi++) {
		  CAPTURE(jphi);
		  CHECK(Approx(q.grad_grad_B_tensor(jphi, j1, j2, j3)).epsilon(tol) == tensor(jphi, j1, j2, j3));

		  // For all configs, the tensor should be symmetric in the 1st 2 indices:
		  CHECK(Approx(q.grad_grad_B_tensor(jphi, j1, j2, j3)).epsilon(tol) == q.grad_grad_B_tensor(jphi, j2, j1, j3));

		  // For curl-free fields, the tensor should also be symmetric in the last 2 indices:
		  if (jconfig == 0 || jconfig == 1 || jconfig == 3) {
		    CHECK(Approx(q.grad_grad_B_tensor(jphi, j1, j2, j3)).epsilon(tol) == q.grad_grad_B_tensor(jphi, j1, j3, j2));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
