#include <iostream>
#include <vector>
#include <string>
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

TEST_CASE("Compare published configurations to fortran version of QSC") {
  std::vector<std::string> cconfigs = {
    "r1 section 5.1",
    "r1 section 5.2",
    "r1 section 5.3",
    "r2 section 5.1",
    "r2 section 5.2",
    "r2 section 5.3",
    "r2 section 5.4",
    "r2 section 5.5"};
  
  std::vector<std::string> fconfigs = {
    "quasisymmetry_out.LandremanSenguptaPlunk_section5.1_order_r1_finite_r_linear.nc",
    "quasisymmetry_out.LandremanSenguptaPlunk_section5.2_order_r1_finite_r_nonlinear.nc",
    "quasisymmetry_out.LandremanSenguptaPlunk_section5.3_order_r1_finite_r_linear.nc",
    "quasisymmetry_out.LandremanSengupta2019_section5.1.nc",
    "quasisymmetry_out.LandremanSengupta2019_section5.2.nc",
    "quasisymmetry_out.LandremanSengupta2019_section5.3.nc",
    "quasisymmetry_out.LandremanSengupta2019_section5.4.nc",
    "quasisymmetry_out.LandremanSengupta2019_section5.5.nc"};
  
  qscfloat tol, well_tol;
  if (single) {
    tol = 3.0e-3;
    well_tol = 3.0e-4;
  } else {
    tol = 1.0e-10;
    well_tol = 1.0e-12;
  }
  
  Qsc f;
  std::cout << "Hello world" << std::endl;
  for (int jconfig = 0; jconfig < fconfigs.size(); jconfig++) {
    std::cout << "jconfig=" << jconfig << std::endl;
    CAPTURE(jconfig);
    f.read_netcdf(fconfigs[jconfig], 'F');
    std::cout << "Read NetCDF fortran file" << std::endl;
    Qsc c(cconfigs[jconfig]);
    c.nphi = f.nphi;
    c.calculate();
    
    // Scalars
    CHECK(c.nphi == f.nphi);
    CHECK(c.nfp == f.nfp);
    CHECK(Approx(c.eta_bar) == f.eta_bar);
    CHECK(Approx(c.B0) == f.B0);
    CHECK(c.sG == f.sG);
    CHECK(c.spsi == f.spsi);
    CHECK(Approx(c.axis_length) == f.axis_length);
    CHECK(Approx(c.abs_G0_over_B0) == f.abs_G0_over_B0);
    CHECK(Approx(c.rms_curvature) == f.rms_curvature);
    CHECK(Approx(c.mean_elongation) == f.mean_elongation);
    CHECK(Approx(c.standard_deviation_of_R) == f.standard_deviation_of_R);
    CHECK(Approx(c.standard_deviation_of_Z) == f.standard_deviation_of_Z);
    CHECK(Approx(c.iota) == f.iota);
    CHECK(Approx(c.sigma0) == f.sigma0);
    CHECK(c.helicity == f.helicity);
    if (c.at_least_order_r2) {
      CHECK(Approx(c.p2) == f.p2);
      CHECK(Approx(c.B2s) == f.B2s);
      CHECK(Approx(c.B2c) == f.B2c);
      CHECK(Approx(c.B20_mean) == f.B20_mean);
      CHECK(Approx(c.B20_residual) == f.B20_residual);
      std::cout << "c.B20_mean:" << c.B20_mean << " f.B20_mean:" << f.B20_mean << " difference:" << c.B20_mean - f.B20_mean << std::endl;
      std::cout << "c.G0:" << c.G0 << " f.abs_G0_over_B0:" << f.abs_G0_over_B0 << " difference:" << c.G0 - f.abs_G0_over_B0 << std::endl;
      CHECK(Approx(c.d2_volume_d_psi2).epsilon(well_tol) == f.d2_volume_d_psi2);
      CHECK(Approx(c.DWell_times_r2) == f.DWell_times_r2);
      CHECK(Approx(c.DGeod_times_r2) == f.DGeod_times_r2);
      CHECK(Approx(c.DMerc_times_r2) == f.DMerc_times_r2);
    }
    // CHECK(Approx(c.) == f.);
    
    // Vectors
    for (int j = 0; j < f.nphi; j++) {
      CHECK(Approx(c.phi[j]) == f.phi[j]);
      CHECK(Approx(c.curvature[j]) == f.curvature[j]);
      CHECK(Approx(c.torsion[j]) == f.torsion[j]);
      CHECK(Approx(c.sigma[j]) == f.sigma[j]);
      CHECK(Approx(c.X1c[j]) == f.X1c[j]);
      CHECK(Approx(c.Y1c[j]) == f.Y1c[j]);
      CHECK(Approx(c.Y1s[j]) == f.Y1s[j]);
      CHECK(Approx(c.R0[j]) == f.R0[j]);
      CHECK(Approx(c.Z0[j]) == f.Z0[j]);
      CHECK(Approx(c.d_l_d_phi[j]) == f.d_l_d_phi[j]);
      CHECK(Approx(c.elongation[j]) == f.elongation[j]);
      CHECK(Approx(c.Boozer_toroidal_angle[j]) == f.Boozer_toroidal_angle[j]);
      CHECK(Approx(c.L_grad_B_inverse[j]) == f.L_grad_B_inverse[j]);
      if (c.at_least_order_r2) {
	// For O(r^2) quantities, we need a loose tolerance for single precision
	CHECK(Approx(c.B20[j]).epsilon(tol) == f.B20[j]);
	CHECK(Approx(c.X20[j]).epsilon(tol) == f.X20[j]);
	CHECK(Approx(c.X2s[j]).epsilon(tol) == f.X2s[j]);
	CHECK(Approx(c.X2c[j]).epsilon(tol) == f.X2c[j]);
	CHECK(Approx(c.Y20[j]).epsilon(tol) == f.Y20[j]);
	CHECK(Approx(c.Y2s[j]).epsilon(tol) == f.Y2s[j]);
	CHECK(Approx(c.Y2c[j]).epsilon(tol) == f.Y2c[j]);
	CHECK(Approx(c.Z20[j]).epsilon(tol) == f.Z20[j]);
	CHECK(Approx(c.Z2s[j]).epsilon(tol) == f.Z2s[j]);
	CHECK(Approx(c.Z2c[j]).epsilon(tol) == f.Z2c[j]);
      }
    }
  }
}
