#include <chrono>
#include <vector>
#include <mpi.h>
#include "multiopt_scan.hpp"
#include "netcdf_writer.hpp"

using namespace qsc;

void MultiOptScan::write_netcdf() {
  // Only proc 0 should run this subroutine:
  if (!proc0) return;
  
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) start = std::chrono::steady_clock::now();

  if (verbose > 0) std::cout << "Writing output to " << outfilename << std::endl;
  qsc::NetCDFWriter nc(outfilename, false);

  // Define dimensions
  dim_id_type nphi_dim, axis_nmax_plus_1_dim, n_scan_dim;
  // It is tricky to write data that depends on the output of a qsc
  //solve, because proc 0 never actually does any solves in a MultiOpt
  //scan.
  // nphi_dim = nc.dim("nphi", mo_ref.opts[0].q.nphi);
  axis_nmax_plus_1_dim = nc.dim("axis_nmax_plus_1", axis_nmax_plus_1);
  n_scan_dim = nc.dim("n_scan", n_scan);
  
  // Scalars
  std::string general_option = "multiopt_scan";
  nc.put("general_option", general_option, "Whether this job was a single configuration vs a scan");
  
  nc.put("nfp", mo_ref.opts[0].q.nfp, "Number of field periods", "dimensionless");
  //nc.put("nphi", mo.opts[0].q.nphi, "Number of grid points in the toroidal angle phi", "dimensionless");
  // In the next line, we cast n_scan to an int because long long ints require netcdf-4, which cannot be read by scipy.io.netcdf.
  
  int n_scan_int = (int)n_scan;
  nc.put("n_scan", n_scan_int, "Number of configurations kept from the scan and saved in this file", "dimensionless");
  nc.put("attempts", filters[ATTEMPTS], "Number of configurations examined in the scan", "dimensionless");
  nc.put("rejected_due_to_R0", filters[REJECTED_DUE_TO_R0], "Number of configurations in the scan that were rejected due to R0 becoming <= 0", "dimensionless");
  nc.put("rejected_due_to_iota", filters[REJECTED_DUE_TO_IOTA], "Number of configurations in the scan that were rejected due to the min_iota_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_elongation", filters[REJECTED_DUE_TO_ELONGATION], "Number of configurations in the scan that were rejected due to the max_elongation_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_L_grad_B", filters[REJECTED_DUE_TO_L_GRAD_B], "Number of configurations in the scan that were rejected due to the min_L_grad_B_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_B20_variation", filters[REJECTED_DUE_TO_B20_VARIATION], "Number of configurations in the scan that were rejected due to the max_B20_variation_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_L_grad_grad_B", filters[REJECTED_DUE_TO_L_GRAD_GRAD_B], "Number of configurations in the scan that were rejected due to the min_L_grad_grad_B_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_d2_volume_d_psi2", filters[REJECTED_DUE_TO_D2_VOLUME_D_PSI2], "Number of configurations in the scan that were rejected due to the max_d2_volume_d_psi2_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_DMerc", filters[REJECTED_DUE_TO_DMERC], "Number of configurations in the scan that were rejected due to the min_DMerc_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_r_singularity", filters[REJECTED_DUE_TO_R_SINGULARITY], "Number of configurations in the scan that were rejected due to the min_r_singularity_to_keep filter", "dimensionless");

  nc.put("fraction_kept", filter_fractions[KEPT], "Fraction of the attempted configurations from the scan that were kept and saved in this file", "dimensionless");
  nc.put("fraction_rejected_due_to_R0", filter_fractions[REJECTED_DUE_TO_R0], "Fraction of configurations in the scan that were rejected due to R0 becoming <= 0", "dimensionless");
  nc.put("fraction_rejected_due_to_iota", filter_fractions[REJECTED_DUE_TO_IOTA], "Fraction of configurations in the scan that were rejected due to the min_iota_to_keep filter", "dimensionless");
  nc.put("fraction_rejected_due_to_elongation", filter_fractions[REJECTED_DUE_TO_ELONGATION], "Fraction of configurations in the scan that were rejected due to the max_elongation_to_keep filter", "dimensionless");
  nc.put("fraction_rejected_due_to_L_grad_B", filter_fractions[REJECTED_DUE_TO_L_GRAD_B], "Fraction of configurations in the scan that were rejected due to the min_L_grad_B_to_keep filter", "dimensionless");
  nc.put("fraction_rejected_due_to_B20_variation", filter_fractions[REJECTED_DUE_TO_B20_VARIATION], "Fraction of configurations in the scan that were rejected due to the max_B20_variation_to_keep filter", "dimensionless");
  nc.put("fraction_rejected_due_to_L_grad_grad_B", filter_fractions[REJECTED_DUE_TO_L_GRAD_GRAD_B], "Fraction of configurations in the scan that were rejected due to the min_L_grad_grad_B_to_keep filter", "dimensionless");
  nc.put("fraction_rejected_due_to_d2_volume_d_psi2", filter_fractions[REJECTED_DUE_TO_D2_VOLUME_D_PSI2], "Fraction of configurations in the scan that were rejected due to the max_d2_volume_d_psi2_to_keep filter", "dimensionless");
  nc.put("fraction_rejected_due_to_DMerc", filter_fractions[REJECTED_DUE_TO_DMERC], "Fraction of configurations in the scan that were rejected due to the min_DMerc_to_keep filter", "dimensionless");
  nc.put("fraction_rejected_due_to_r_singularity", filter_fractions[REJECTED_DUE_TO_R_SINGULARITY], "Fraction of configurations in the scan that were rejected due to the min_r_singularity_to_keep filter", "dimensionless");

  int keep_all_int = (int) keep_all;
  nc.put("keep_all", keep_all_int, "1 if all configurations from the scan were saved, 0 if some configurations were filtered out", "dimensionless");
;
  if (!keep_all) {
    nc.put("min_R0_to_keep", min_R0_to_keep, "Configurations were kept in the scan only if the major radius of the magnetic axis was at least this value", "meter");
    nc.put("min_iota_to_keep", min_iota_to_keep, "Configurations were kept in the scan only if the absolute value of the on-axis rotational transform was at least this value", "dimensionless");
    nc.put("max_elongation_to_keep", max_elongation_to_keep, "Configurations were kept in the scan only if the elongation (in the plane perpendicular to the magnetic axis) was no greater than this value at all toroidal angles", "dimensionless");
    nc.put("min_L_grad_B_to_keep", min_L_grad_B_to_keep, "Configurations were kept in the scan only if the scale length L_grad_B (eq (3.1) in Landreman J Plasma Physics (2021)) is at least this value at each toroidal angle", "meter");
    nc.put("min_L_grad_grad_B_to_keep", min_L_grad_grad_B_to_keep, "Configurations were kept in the scan only if the scale length L_grad_grad_B (eq (3.2) in Landreman J Plasma Physics (2021)) is at least this value at each toroidal angle", "meter");
    nc.put("max_B20_variation_to_keep", max_B20_variation_to_keep, "Configurations were kept in the scan only if the range (maximum minus minimum) of B20 over toroidal angle is no more than this value", "Tesla/(meter^2)");
    nc.put("min_r_singularity_to_keep", min_r_singularity_to_keep, "Configurations were kept in the scan only if r_singularity_robust is at least this value. r_singularity_robust is the robust estimate of the minor radius at which the flux surface shapes become singular, r_c, as detailed in section 4.2 of Landreman, J Plasma Physics (2021)", "meter");
    nc.put("max_d2_volume_d_psi2_to_keep", max_d2_volume_d_psi2_to_keep, "Configurations were kept in the scan only if the magnetic well d2_volume_d_psi2 is no more than this value. d2_volume_d_psi2 is the second derivative of flux surface volume with respect to psi, where 2*pi*psi is the toroidal flux.", "Tesla^{-2} meter^{-1}");
    nc.put("min_DMerc_to_keep", min_DMerc_to_keep, "Configurations were kept in the scan only if the Mercier stability criterion DMerc_times_r2 is at least this value. DMerc_times_r2 corresponds to r^2 times the quantity DMerc in Landreman and Jorge, J Plasma Phys (2020).", "Tesla^{-2} meter^{-2}");
  }
  nc.put("max_newton_iterations", mo_ref.opts[0].q.max_newton_iterations, "Maximum iterations of Newton's method for solving the sigma equation", "dimensionless");
  nc.put("max_linesearch_iterations", mo_ref.opts[0].q.max_linesearch_iterations, "Maximum number of times the step size is reduced in the line search for each iteration of Newton's method when solving the sigma equation", "dimensionless");
  nc.put("newton_tolerance", mo_ref.opts[0].q.newton_tolerance, "L2 norm of the residual used as a stopping criterion for Newton's method when solving the sigma equation", "dimensionless");
  nc.put("I2", mo_ref.opts[0].q.I2, "r^2 term in I(r), which is the toroidal current inside the flux surface times mu0/(2pi)", "Tesla/meter");
  //nc.put("d_phi", mo_ref.opts[0].q.d_phi, "Grid spacing in phi", "dimensionless");
  nc.put("B0", mo_ref.opts[0].q.B0, "Magnetic field magnitude on the magnetic axis", "Telsa");
  nc.put("sG", mo_ref.opts[0].q.sG, "Sign of G0", "dimensionless");
  nc.put("spsi", mo_ref.opts[0].q.spsi, "Sign of the toroidal flux psi", "dimensionless");
  nc.put("p2", mo_ref.opts[0].q.p2, "r^2 term in p(r), the pressure profile", "Pascal/(meter^2)");

  // 1D arrays
  //nc.put(nphi_dim, "phi", mo_ref.opts[0].q.phi, "The grid in the standard toroidal angle phi", "dimensionless");
  nc.put(n_scan_dim, "scan_eta_bar", scan_eta_bar, "For each configuration kept from the scan, the constant equal to B1c / B0", "1/meter");
  nc.put(n_scan_dim, "scan_sigma0", scan_sigma0, "For each configuration kept from the scan, the value of sigma at phi=0", "dimensionless");
  nc.put(n_scan_dim, "scan_B2c", scan_B2c, "For each configuration kept from the scan, the r^2 * cos(2*theta) term in |B|", "Tesla/(meter^2)");
  nc.put(n_scan_dim, "scan_B2s", scan_B2s, "For each configuration kept from the scan, the r^2 * sin(2*theta) term in |B|", "Tesla/(meter^2)");
  nc.put(n_scan_dim, "scan_iota", scan_iota, "For each configuration kept from the scan, the rotational transform on axis", "dimensionless");
  nc.put(n_scan_dim, "scan_min_R0", scan_min_R0, "For each configuration kept from the scan, the minimum value of R0, the major radius of the magnetic axis. This variable corresponds to grid_min_R0 in a single Qsc calculation.", "meter");
  nc.put(n_scan_dim, "scan_max_curvature", scan_max_curvature, "For each configuration kept from the scan, the maximum curvature of the magnetic axis", "1/meter");
  nc.put(n_scan_dim, "scan_max_elongation", scan_max_elongation, "For each configuration kept from the scan, the maximum along the magnetic axis of the elongation in the plane perpendicular to the axis", "dimensionless");
  nc.put(n_scan_dim, "scan_min_L_grad_B", scan_min_L_grad_B, "For each configuration kept from the scan, the minimum along the magnetic axis of the scale length L_grad_B, (eq (3.1) in Landreman J Plasma Physics (2021). This quantity corresponds to grid_min_L_grad_B for a single Qsc run.", "meter");
  nc.put(n_scan_dim, "scan_helicity", scan_helicity, "For each configuration kept from the scan, the number of times the normal vector of the magnetic axis rotates poloidally as the axis is followed toroidally for one field period. The integer N appearing in our papers is equal to -helicity * nfp.", "dimensionless");
  nc.put(n_scan_dim, "scan_standard_deviation_of_R", scan_standard_deviation_of_R, "Standard deviation of the major radius of the magnetic axis, with respect to arclength along the axis", "meter");
  nc.put(n_scan_dim, "scan_standard_deviation_of_Z", scan_standard_deviation_of_Z, "Standard deviation of the Cartesian Z coordinate of the magnetic axis, with respect to arclength along the axis", "meter");
  nc.put(n_scan_dim, "scan_min_L_grad_grad_B", scan_min_L_grad_grad_B, "For each configuration kept from the scan, the minimum along the magnetic axis of the scale length L_grad_grad_B, (eq (3.2) in Landreman J Plasma Physics (2021). This quantity corresponds to grid_min_L_grad_grad_B for a single Qsc run.", "meter");
  nc.put(n_scan_dim, "scan_B20_variation", scan_B20_variation, "For each configuration kept from the scan, the maximum of B20 along the magnetic axis minus the minimum of B20. This quantity corresponds to B20_grid_variation for a single Qsc run.", "Telsa/(meter^2)");
  nc.put(n_scan_dim, "scan_B20_residual", scan_B20_residual, "", "");
  nc.put(n_scan_dim, "scan_r_singularity", scan_r_singularity, "For each configuration kept from the scan, the value of r_singularity_robust. r_singularity_robust is the robust estimate of the minor radius at which the flux surface shapes become singular, r_c, as detailed in section 4.2 of Landreman, J Plasma Physics (2021)", "meter");
  nc.put(n_scan_dim, "scan_d2_volume_d_psi2", scan_d2_volume_d_psi2, "For each configuration kept from the scan, the value of magnetic well d2_volume_d_psi2, the second derivative of flux surface volume with respect to psi, where 2*pi*psi is the toroidal flux.", "Tesla^{-2} meter^{-1}");
  nc.put(n_scan_dim, "scan_DMerc_times_r2", scan_DMerc_times_r2, "For each configuration kept from the scan, the overall Mercier stability criterion times the square of the effective minor radius r. This quantity corresponds to DMerc_times_r2 for a single Qsc run. DMerc (without the r^2) corresponds to the quantity DMerc in VMEC, and to DMerc in Landreman and Jorge, J Plasma Phys (2020).", "Tesla^{-2} meter^{-2}");

  nc.put(n_scan_dim, "scan_weight_B20", scan_weight_B20, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_iota", scan_weight_iota, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_target_iota", scan_target_iota, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_elongation", scan_weight_elongation, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_curvature", scan_weight_curvature, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_R0", scan_weight_R0, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_target_min_R0", scan_target_min_R0, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_d2_volume_d_psi2", scan_weight_d2_volume_d_psi2, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_max_d2_volume_d_psi2", scan_max_d2_volume_d_psi2, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_XY2", scan_weight_XY2, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_XY2Prime", scan_weight_XY2Prime, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_Z2", scan_weight_Z2, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_Z2Prime", scan_weight_Z2Prime, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_XY3", scan_weight_XY3, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_XY3Prime", scan_weight_XY3Prime, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_grad_B", scan_weight_grad_B, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_grad_grad_B", scan_weight_grad_grad_B, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_r_singularity", scan_weight_r_singularity, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_axis_length", scan_weight_axis_length, " ", "dimensionless");
  nc.put(n_scan_dim, "scan_weight_standard_deviation_of_R", scan_weight_standard_deviation_of_R, " ", "dimensionless");

  // ND arrays for N > 1:
  std::vector<dim_id_type> axis_nmax_plus_1_n_scan_dim {axis_nmax_plus_1_dim, n_scan_dim};
  nc.put(axis_nmax_plus_1_n_scan_dim, "scan_R0c", &scan_R0c(0, 0), "For each configuration kept from the scan, the amplitudes of the cos(n*phi) components of the major radius of the magnetic axis", "meter");
  nc.put(axis_nmax_plus_1_n_scan_dim, "scan_R0s", &scan_R0s(0, 0), "For each configuration kept from the scan, the amplitudes of the sin(n*phi) components of the major radius of the magnetic axis", "meter");
  nc.put(axis_nmax_plus_1_n_scan_dim, "scan_Z0c", &scan_Z0c(0, 0), "For each configuration kept from the scan, the amplitudes of the cos(n*phi) components of the Cartesian Z coordinate of the magnetic axis", "meter");
  nc.put(axis_nmax_plus_1_n_scan_dim, "scan_Z0s", &scan_Z0s(0, 0), "For each configuration kept from the scan, the amplitudes of the sin(n*phi) components of the Cartesian Z coordinate of the magnetic axis", "meter");

  // Done defining the NetCDF data.
  nc.write_and_close();
  
  if (verbose > 0) {
    auto end = std::chrono::steady_clock::now();    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for write_netcdf: "
              << elapsed.count() << " seconds" << std::endl;
  }
}
