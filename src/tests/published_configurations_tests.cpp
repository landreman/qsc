#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("Check iota for published configurations") {
  int nphis[] = {50, 63};
  int nphi;
  qscfloat loose = 0.01; // Tolerance for comparing grid_max vs max quantities.
  
  for (int j = 0; j < 2; j++) {
    nphi = nphis[j];
    CAPTURE(nphi);
    
    Qsc q11("r1 section 5.1");
    q11.nphi = nphi;
    q11.calculate();
    CHECK(Approx(q11.iota) == 0.418306910215178);
    CHECK(Approx(q11.mean_elongation) == 2.28434216811829);
    CHECK(Approx(q11.grid_max_elongation).epsilon(loose) == 2.41373705531443);
    CHECK(Approx(q11.grid_min_L_grad_B).epsilon(loose) == 1 / 1.52948586064743);
    
    Qsc q12("r1 section 5.2");
    q12.nphi = nphi;
    q12.calculate();
    CHECK(Approx(q12.iota) == 1.93109725535729);
    CHECK(Approx(q12.mean_elongation) == 2.12218817610318);
    CHECK(Approx(q12.grid_max_elongation).epsilon(loose) == 3.08125973323805);
    CHECK(Approx(q12.grid_min_L_grad_B).epsilon(loose) == 1 / 4.73234243198959);
    
    Qsc q13("r1 section 5.3");
    q13.nphi = nphi;
    q13.calculate();
    CHECK(Approx(q13.iota) == 0.311181373123728);
    CHECK(Approx(q13.mean_elongation) == 2.48657801778199);
    CHECK(Approx(q13.grid_max_elongation).epsilon(loose) == 3.30480616121377);
    CHECK(Approx(q13.grid_min_L_grad_B).epsilon(loose) == 1 / 1.7014044379421);
    
    Qsc q21("r2 section 5.1");
    q21.nphi = nphi;
    q21.calculate();
    CHECK(Approx(q21.iota) == -0.420473351810416);
    CHECK(Approx(q21.mean_elongation) == 3.58268292490318);
    CHECK(Approx(q21.grid_max_elongation).epsilon(loose) == 4.38384260252044);
    CHECK(Approx(q21.grid_min_L_grad_B).epsilon(loose) == 1 / 1.39153088147691);
    
    Qsc q22("r2 section 5.2");
    q22.nphi = nphi;
    q22.calculate();
    CHECK(Approx(q22.iota) == -0.423723995700502);
    CHECK(Approx(q22.mean_elongation) == 3.61629912951486);
    CHECK(Approx(q22.grid_max_elongation).epsilon(loose) == 4.86202324600918);
    CHECK(Approx(q22.grid_min_L_grad_B).epsilon(loose) == 1 / 1.47675199709439);
    
    Qsc q23("r2 section 5.3");
    q23.nphi = nphi;
    q23.calculate();
    CHECK(Approx(q23.iota) == 0.959698159859113);
    CHECK(Approx(q23.mean_elongation) == 1.8447534972894);
    CHECK(Approx(q23.grid_max_elongation).epsilon(loose) == 2.20914173760329);
    CHECK(Approx(q23.grid_min_L_grad_B).epsilon(loose) == 1 / 1.4922510395338);
    
    Qsc q24("r2 section 5.4");
    q24.nphi = nphi;
    q24.calculate();
    CHECK(Approx(q24.iota) == -1.14413695118515);
    CHECK(Approx(q24.mean_elongation) == 2.87255662325544);
    CHECK(Approx(q24.grid_max_elongation).epsilon(loose) == 2.98649978627541);
    CHECK(Approx(q24.grid_min_L_grad_B).epsilon(loose) == 1 / 2.64098280647292);
    
    Qsc q25("r2 section 5.5");
    q25.nphi = nphi;
    q25.calculate();
    CHECK(Approx(q25.iota) == -0.828885267089981);
    CHECK(Approx(q25.mean_elongation) == 2.14382600115829);
    CHECK(Approx(q25.grid_max_elongation).epsilon(loose) == 3.6226360623368);
    CHECK(Approx(q25.grid_min_L_grad_B).epsilon(loose) == 1 / 4.85287603883526);
  }
}

