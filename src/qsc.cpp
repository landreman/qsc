#include <ctime>
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

Qsc::Qsc() :
  // Call constructor of member objects:
  d_d_phi(1, 1),
  d_d_varphi(1, 1),
  work_matrix(1, 1),
  grad_B_tensor(1, 3, 3),
  grad_grad_B_tensor(1, 3, 3, 3),
  r2_matrix(1, 1)
{
  defaults();
  
  R0c.resize(1, 1.0);
  R0s.resize(1, 0.0);
  Z0c.resize(1, 0.0);
  Z0s.resize(1, 0.0);

}

/** High-level routine to call the low-level routines.
 */
void Qsc::calculate() {
  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }

  validate();
  allocate();
  init_axis();
  solve_sigma_equation();
  r1_diagnostics();
  if (at_least_order_r2) {
    calculate_r2();
    r2_diagnostics();
  }

  if (verbose > 0) {
    end_time = std::clock();
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for calculate() from chrono:           "
              << elapsed.count() << " seconds" << std::endl;
    std::cout << "Time for calculate() from ctime (CPU time): "
              << double(end_time - start_time) / CLOCKS_PER_SEC
              << " seconds" << std::endl;
  }

}

void Qsc::run(std::string directory_and_infile) {
  std::string directory_and_outfile = qsc::outfile(directory_and_infile);
  input(directory_and_infile);
  calculate();
  write_netcdf(directory_and_outfile);
}
