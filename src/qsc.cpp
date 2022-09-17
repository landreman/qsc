#include <chrono>
#include "qsc.hpp"

using namespace qsc;

void Qsc::defaults() {
  // Set defaults.
  verbose = 1;
  order_r_option = ORDER_R_OPTION_R1;
  
  sG = 1;
  spsi = 1;
  B0 = 1.0;
  eta_bar = -1.0;
  I2 = 0.0;
  sigma0 = 0.0;
  B2s = 0.0;
  B2c = 0.0;
  p2 = 0.0;

  nfp = 3;
  nphi = 15;

  max_newton_iterations = 12;
  max_linesearch_iterations = 4;
  if (single) {
    newton_tolerance = 1.0e-5;
  } else {
    newton_tolerance = 1.0e-12;
  }

  order_r_option = "r1";
}

Qsc::Qsc() {
  defaults();

  resize_axis_arrays(1, 0.0);
  R0c[0] = 1.0;

  grid_min_R0 = 0.0;
  curvature = 0.0;
  grid_max_curvature = 0.0;
  iota = 0.0;
  elongation = 0.0;
  grid_max_elongation = 0.0;
  L_grad_B = 0.0;
  grid_min_L_grad_B = 0.0;
  L_grad_grad_B = 0.0;
  grid_min_L_grad_grad_B = 0.0;
  r_singularity_robust = 0.0;
  r_hat_singularity_robust = 0.0;
  helicity = 0;
  B20_grid_variation = 0.0;
  B20_residual = 0.0;
  d2_volume_d_psi2 = 0.0;
  DMerc_times_r2 = 0.0;
  DGeod_times_r2 = 0.0;
  DWell_times_r2 = 0.0;
}

void Qsc::resize_axis_arrays(int newsize, qscfloat val) {
  R0c.resize(newsize, val);
  R0s.resize(newsize, val);
  Z0c.resize(newsize, val);
  Z0s.resize(newsize, val);
  fc.resize(newsize, val);
  fs.resize(newsize, val);
}

/** High-level routine to call the low-level routines.
 */
void Qsc::init() {
  validate();
  allocate();
}

/** High-level routine to call the low-level routines.
 */
void Qsc::calculate() {
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) start = std::chrono::steady_clock::now();

  init_axis();
  solve_sigma_equation();
  r1_diagnostics();
  if (at_least_order_r2) {
    calculate_r2();
    r2_diagnostics();
  }

  if (verbose > 0) {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for calculate(): "
              << elapsed.count() << " seconds" << std::endl;
  }

}

void Qsc::run(std::string directory_and_infile) {
  std::string directory_and_outfile = qsc::outfile(directory_and_infile);
  input(directory_and_infile);
  init();
  calculate();
  write_netcdf(directory_and_outfile);
}
