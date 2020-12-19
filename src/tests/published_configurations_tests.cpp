#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("Check iota for published configurations") {
  int nphis[] = {50, 63};
  int nphi;
  for (int j = 0; j < 2; j++) {
    nphi = nphis[j];
    CAPTURE(nphi);
    
    Qsc q11("r1 section 5.1");
    q11.nphi = nphi;
    q11.calculate();
    CHECK(Approx(q11.iota) == 0.418306910215178);
    CHECK(Approx(q11.mean_elongation) == 2.28434216811829);
    
    Qsc q12("r1 section 5.2");
    q12.nphi = nphi;
    q12.calculate();
    CHECK(Approx(q12.iota) == 1.93109725535729);
    CHECK(Approx(q12.mean_elongation) == 2.12218817610318);
    
    Qsc q13("r1 section 5.3");
    q13.nphi = nphi;
    q13.calculate();
    CHECK(Approx(q13.iota) == 0.311181373123728);
    CHECK(Approx(q13.mean_elongation) == 2.48657801778199);
    
    Qsc q21("r2 section 5.1");
    q21.nphi = nphi;
    q21.calculate();
    CHECK(Approx(q21.iota) == -0.420473351810416);
    CHECK(Approx(q21.mean_elongation) == 3.58268292490318);
    
    Qsc q22("r2 section 5.2");
    q22.nphi = nphi;
    q22.calculate();
    CHECK(Approx(q22.iota) == -0.423723995700502);
    CHECK(Approx(q22.mean_elongation) == 3.61629912951486);
    
    Qsc q23("r2 section 5.3");
    q23.nphi = nphi;
    q23.calculate();
    CHECK(Approx(q23.iota) == 0.959698159859113);
    CHECK(Approx(q23.mean_elongation) == 1.8447534972894);
    
    Qsc q24("r2 section 5.4");
    q24.nphi = nphi;
    q24.calculate();
    CHECK(Approx(q24.iota) == -1.14413695118515);
    CHECK(Approx(q24.mean_elongation) == 2.87255662325544);
    
    Qsc q25("r2 section 5.5");
    q25.nphi = nphi;
    q25.calculate();
    CHECK(Approx(q25.iota) == -0.828885267089981);
    CHECK(Approx(q25.mean_elongation) == 2.14382600115829);
  }
}

