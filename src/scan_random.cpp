#include "qsc.hpp"
#include "scan.hpp"
#include "random.hpp"

using namespace qsc;

void Scan::random() {
  const int n_parameters = 14;
  const int n_int_parameters = 1;
  const int axis_nmax_plus_1 = R0c_max.size();
  Matrix parameters_local(max_keep_per_proc, n_parameters);
  Matrix int_parameters_local(max_keep_per_proc, n_int_parameters);
  Matrix fourier_parameters_local(axis_nmax_plus_1 * 4, n_parameters);
  big j_scan = 0, attempts, attempts_local = 0, n_kept;
  big rejected_due_to_R0_crude = 0, rejected_due_to_R0 = 0, rejected_due_to_curvature = 0;
  big rejected_due_to_iota = 0, rejected_due_to_elongation = 0;
  big rejected_due_to_L_grad_B = 0, rejected_due_to_L_grad_grad_B = 0;
  big rejected_due_to_B20_variation = 0, rejected_due_to_r_singularity = 0;
  big rejected_due_to_d2_volume_d_psi2 = 0, rejected_due_to_DMerc = 0;
  int j;
  qscfloat R0_at_0, R0_at_half_period, val;
  
  // Initialize random distributions:
  Random random_eta_bar(deterministic, eta_bar_scan_option, eta_bar_min, eta_bar_max);
  Random random_sigma0(deterministic, sigma0_scan_option, sigma0_min, sigma0_max);
  Random random_B2c(deterministic, B2c_scan_option, B2c_min, B2c_max);
  Random random_B2s(deterministic, B2s_scan_option, B2s_min, B2s_max);
  Random* random_R0c[axis_nmax_plus_1];
  Random* random_R0s[axis_nmax_plus_1];
  Random* random_Z0c[axis_nmax_plus_1];
  Random* random_Z0s[axis_nmax_plus_1];
  for (j = 0; j < axis_nmax_plus_1; j++) {
    random_R0c[j] = new Random(deterministic, fourier_scan_option, R0c_min[j], R0c_max[j]);
    random_R0s[j] = new Random(deterministic, fourier_scan_option, R0s_min[j], R0s_max[j]);
    random_Z0c[j] = new Random(deterministic, fourier_scan_option, Z0c_min[j], Z0c_max[j]);
    random_Z0s[j] = new Random(deterministic, fourier_scan_option, Z0s_min[j], Z0s_max[j]);
  }
  
  std::chrono::time_point<std::chrono::steady_clock> start_time, end_time;
  start_time = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed;

  q.validate();
  q.allocate();
  
  bool keep_going = true;
  while (keep_going) {
    end_time = std::chrono::steady_clock::now();
    elapsed = end_time - start_time;
    if (elapsed.count() > max_seconds) keep_going = false;

    attempts_local++;

    // Initialize R0c, and do a crude check of whether R0 goes negative:
    R0_at_half_period = 0;
    for (j = 0; j < axis_nmax_plus_1; j++) {
      val = random_R0c[j]->get();
      q.R0c[j] = val;
      if (j % 2 == 0) {
	R0_at_half_period += val;
      } else {
	R0_at_half_period -= val;
      }
    }
    R0_at_0 = q.R0c.sum();
    if (R0_at_0 <= 0 || R0_at_half_period <= 0) {
      rejected_due_to_R0_crude++;
      continue;
    }

    // Initialize the rest of the axis shape:
    for (j = 0; j < axis_nmax_plus_1; j++) {
      q.R0s[j] = random_R0s[j]->get();
      q.Z0s[j] = random_Z0s[j]->get();
      q.Z0c[j] = random_Z0s[j]->get();
    }
    q.init_axis();
    if (q.grid_min_R0 < min_R0_to_keep) {
      rejected_due_to_R0++;
      continue;
    }
    if (!keep_all && 1.0 / q.grid_max_curvature < min_L_grad_B_to_keep) {
      rejected_due_to_curvature++;
      continue;
    }

    // Initialize the rest of the O(r^1) parameters:
    q.eta_bar = random_eta_bar.get();
    q.sigma0 = random_sigma0.get();

    // Here is the main O(r^1) solve:
    q.solve_sigma_equation();
    q.r1_diagnostics();
    if (!keep_all && q.grid_max_elongation > max_elongation_to_keep) {
      rejected_due_to_elongation++;
      continue;
    }
    if (!keep_all && q.grid_min_L_grad_B < min_L_grad_B_to_keep) {
      rejected_due_to_L_grad_B++;
      continue;
    }

    if (q.at_least_order_r2) {
      // Initialize the rest of the O(r^2) parameters:
      q.B2c = random_B2c.get();
      q.B2s = random_B2s.get();

      // Here is the main O(r^2) solve:
      q.calculate_r2();

      // Filter results:
      if (!keep_all && q.B20_grid_variation > max_B20_variation_to_keep) {
	rejected_due_to_B20_variation++;
	continue;
      }

      q.mercier();
      if (!keep_all && q.d2_volume_d_psi2 > max_d2_volume_d_psi2_to_keep) {
	rejected_due_to_d2_volume_d_psi2++;
	continue;
      }
      if (!keep_all && q.DMerc_times_r2 < min_DMerc_to_keep) {
	rejected_due_to_DMerc++;
	continue;
      }

      q.calculate_grad_grad_B_tensor();
      if (!keep_all && q.grid_min_L_grad_grad_B < min_L_grad_grad_B_to_keep) {
	rejected_due_to_L_grad_grad_B++;
	continue;
      }

      q.calculate_r_singularity();
      if (!keep_all && q.r_singularity_robust < min_r_singularity_to_keep) {
	rejected_due_to_r_singularity++;
	continue;
      }
    } // if at_least_order_r2

    // If we made it this far, then we found a keeper.
    parameters_local(j_scan, 0 ) = q.eta_bar;
    parameters_local(j_scan, 1 ) = q.sigma0;
    parameters_local(j_scan, 2 ) = q.B2c;
    parameters_local(j_scan, 3 ) = q.B2s;
    parameters_local(j_scan, 4 ) = q.grid_min_R0;
    parameters_local(j_scan, 5 ) = q.grid_max_curvature;
    parameters_local(j_scan, 6 ) = q.iota;
    parameters_local(j_scan, 7 ) = q.grid_max_elongation;
    parameters_local(j_scan, 8 ) = q.grid_min_L_grad_B;
    parameters_local(j_scan, 9 ) = q.grid_min_L_grad_grad_B;
    parameters_local(j_scan, 10) = q.r_singularity_robust;
    parameters_local(j_scan, 11) = q.d2_volume_d_psi2;
    parameters_local(j_scan, 12) = q.DMerc_times_r2;
    parameters_local(j_scan, 13) = q.B20_grid_variation;

    int_parameters_local(j_scan, 0) = q.helicity;
    
    j_scan++;

    if (j_scan >= max_keep_per_proc) keep_going = false;
  }

  big total_rejected = rejected_due_to_R0_crude
    + rejected_due_to_R0
    + rejected_due_to_curvature
    + rejected_due_to_iota
    + rejected_due_to_elongation
    + rejected_due_to_L_grad_B
    + rejected_due_to_B20_variation
    + rejected_due_to_L_grad_grad_B
    + rejected_due_to_d2_volume_d_psi2
    + rejected_due_to_DMerc
    + rejected_due_to_r_singularity;
  
  std::cout << "Summary of scan results:" << std::endl;
  std::cout << "  Configurations attempted:          " << attempts_local << std::endl;
  std::cout << "  Rejected due to crude R0 check:    " << rejected_due_to_R0_crude << std::endl;
  std::cout << "  Rejected due to min_R0:            " << rejected_due_to_R0 << std::endl;
  std::cout << "  Rejected due to max curvature:     " << rejected_due_to_curvature << std::endl;
  std::cout << "  Rejected due to min iota:          " << rejected_due_to_iota << std::endl;
  std::cout << "  Rejected due to max elongation:    " << rejected_due_to_elongation << std::endl;
  std::cout << "  Rejected due to min L_grad_B:      " << rejected_due_to_L_grad_B << std::endl;
  std::cout << "  Rejected due to B20 variation:     " << rejected_due_to_B20_variation << std::endl;
  std::cout << "  Rejected due to min L_grad_grad_B: " << rejected_due_to_L_grad_grad_B << std::endl;
  std::cout << "  Rejected due to d2_volume_d_psi2:  " << rejected_due_to_d2_volume_d_psi2 << std::endl;
  std::cout << "  Rejected due to DMerc:             " << rejected_due_to_DMerc << std::endl;
  std::cout << "  Rejected due to r_singularity:     " << rejected_due_to_r_singularity << std::endl;
  std::cout << "  Kept:                              " << j_scan << std::endl;
  std::cout << "  Kept + rejected:                   " << j_scan + total_rejected << std::endl;
}
