#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

/** Test that the curvature and torsion match an independent
 *  calculation using the fortran code.
*/
TEST_CASE("curvature and torsion: stellarator-symmetric case") {
  Qsc q;
  q.nfp = 5;
  q.nphi = 15;
  
  q.R0c.resize(4, 0.0);
  q.R0s.resize(4, 0.0);
  q.Z0c.resize(4, 0.0);
  q.Z0s.resize(4, 0.0);

  q.R0c[0] = 1.3;
  q.R0c[1] = 0.3;
  q.R0c[2] = 0.01;
  q.R0c[3] = -0.001;

  q.Z0s[1] = 0.4;
  q.Z0s[2] = -0.02;
  q.Z0s[3] = -0.003;

  q.allocate();
  q.init_axis();

  qscfloat curvature_fortran[] = {1.74354628565018, 1.61776632275718, 1.5167042487094,
    1.9179603622369, 2.95373444883134, 3.01448808361584, 1.7714523990583,
    1.02055493647363, 1.02055493647363, 1.77145239905828, 3.01448808361582,
    2.95373444883135, 1.91796036223691, 1.5167042487094, 1.61776632275717};

  qscfloat torsion_fortran[] = {0.257226801231061, -0.131225053326418, -1.12989287766591,
    -1.72727988032403, -1.48973327005739, -1.34398161921833,
    -1.76040161697108, -2.96573007082039, -2.96573007082041,
    -1.7604016169711, -1.34398161921833, -1.48973327005739,
    -1.72727988032403, -1.12989287766593, -0.13122505332643};

  qscfloat varphi_fortran[] = {0, 0.0909479184372571, 0.181828299105257,
    0.268782689120682, 0.347551637441381, 0.42101745128188,
    0.498195826255542, 0.583626271820683, 0.673010789615233,
    0.758441235180374, 0.835619610154036, 0.909085423994535,
    0.987854372315234, 1.07480876233066, 1.16568914299866};
    
  for (int j = 0; j < q.nphi; j++) {
    CHECK(q.curvature[j] == Approx(curvature_fortran[j]));
    CHECK(q.torsion[j] == Approx(torsion_fortran[j]));
    CHECK(q.Boozer_toroidal_angle[j] == Approx(varphi_fortran[j]));
  }
}

TEST_CASE("curvature and torsion: non-stellarator-symmetric case") {
  Qsc q;
  q.nfp = 5;
  q.nphi = 15;
  
  q.R0c.resize(4, 0.0);
  q.R0s.resize(4, 0.0);
  q.Z0c.resize(4, 0.0);
  q.Z0s.resize(4, 0.0);

  q.R0c[0] = 1.3;
  q.R0c[1] = 0.3;
  q.R0c[2] = 0.01;
  q.R0c[3] = -0.001;

  q.Z0s[1] = 0.4;
  q.Z0s[2] = -0.02;
  q.Z0s[3] = -0.003;

  q.R0s[1] = -0.1;
  q.R0s[2] = -0.03;
  q.R0s[3] = 0.002;

  q.Z0c[0] = 0.3;
  q.Z0c[1] = 0.2;
  q.Z0c[2] = 0.04;
  q.Z0c[3] = 0.004;

  q.allocate();
  q.init_axis();

  qscfloat curvature_fortran[] = {2.10743037699653, 2.33190181686696, 1.83273654023051,
    1.81062232906827, 2.28640008392347, 1.76919841474321, 0.919988560478029,
    0.741327470169023, 1.37147330126897, 2.64680884158075, 3.39786486424852,
    2.47005615416209, 1.50865425515356, 1.18136509189105, 1.42042418970102};

  qscfloat torsion_fortran[] = {-0.167822738386845, -0.0785778346620885, -1.02205137493593,
    -2.05213528002946, -0.964613202459108, -0.593496282035916,
    -2.15852857178204, -3.72911055219339, -1.9330792779459,
    -1.53882290974916, -1.42156496444929, -1.11381642382793,
    -0.92608309386204, -0.868339812017432, -0.57696266498748};

  qscfloat varphi_fortran[] = {0, 0.084185130335249, 0.160931495903817,
    0.232881563535092, 0.300551168190665, 0.368933497012765,
    0.444686439112853, 0.528001290336008, 0.612254611059372,
    0.691096975269652, 0.765820243301147, 0.846373713025902,
    0.941973362938683, 1.05053459351092, 1.15941650366667};
  
  for (int j = 0; j < q.nphi; j++) {
    CHECK(q.curvature[j] == Approx(curvature_fortran[j]));
    CHECK(q.torsion[j] == Approx(torsion_fortran[j]));
    CHECK(q.Boozer_toroidal_angle[j] == Approx(varphi_fortran[j]));
  }
}
