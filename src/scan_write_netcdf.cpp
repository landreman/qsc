#include <ctime>
#include <chrono>
#include <vector>
#include <mpi.h>
#include "qsc.hpp"
#include "scan.hpp"
#include "netcdf_writer.hpp"

using namespace qsc;

void Scan::write_netcdf() {
  // Only proc 0 should run this subroutine:
  int mpi_rank;
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  if (mpi_rank != 0) return;
  
  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }

  if (verbose > 0) std::cout << "Writing output to " << outfilename << std::endl;
  qsc::NetCDFWriter nc(outfilename);

  // Define dimensions
  dim_id_type nphi_dim, axis_nmax_plus_1_dim, n_scan_dim;
  nphi_dim = nc.dim("nphi", q.nphi);
  axis_nmax_plus_1_dim = nc.dim("axis_nmax_plus_1", R0c_max.size());
  n_scan_dim = nc.dim("n_scan", n_scan);
  
  // Scalars
  int at_least_order_r2_int = (int) q.at_least_order_r2;
  nc.put("at_least_order_r2", at_least_order_r2_int, "1 if the O(r^2) equations were solved, 0 if not", "dimensionless");
  nc.put("nfp", q.nfp, "Number of field periods", "dimensionless");
  nc.put("nphi", q.nphi, "Number of grid points in the toroidal angle phi", "dimensionless");
  // In the next line, we cast n_scan to an int because long long ints require netcdf-4, which cannot be read by scipy.io.netcdf.
  int n_scan_int = (int)n_scan;
  nc.put("n_scan", n_scan_int, "Number of configurations kept from the scan and saved in this file", "dimensionless");
  nc.put("attempts", filters[ATTEMPTS], "Number of configurations examined in the scan", "dimensionless");
  nc.put("n_sigma_eq_solves", filters[N_SIGMA_EQ_SOLVES], "Number of times the sigma equation was solved during the scan", "dimensionless");
  nc.put("n_r2_solves", filters[N_R2_SOLVES], "Number of times the O(r^2) equations were solved during the scan", "dimensionless");
  nc.put("rejected_due_to_R0_crude", filters[REJECTED_DUE_TO_R0_CRUDE], "Number of configurations in the scan that were rejected due to R0 becoming <= 0 at toroidal angle 0 or half field period", "dimensionless");
  nc.put("rejected_due_to_R0", filters[REJECTED_DUE_TO_R0], "Number of configurations in the scan that were rejected due to R0 becoming <= 0", "dimensionless");
  nc.put("rejected_due_to_curvature", filters[REJECTED_DUE_TO_CURVATURE], "Number of configurations in the scan that were rejected due to the curvature of the magnetic axis exceeding 1 / min_L_grad_B_to_keep", "dimensionless");
  nc.put("rejected_due_to_iota", filters[REJECTED_DUE_TO_IOTA], "Number of configurations in the scan that were rejected due to the min_iota_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_elongation", filters[REJECTED_DUE_TO_ELONGATION], "Number of configurations in the scan that were rejected due to the max_elongation_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_L_grad_B", filters[REJECTED_DUE_TO_L_GRAD_B], "Number of configurations in the scan that were rejected due to the min_L_grad_B_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_B20_variation", filters[REJECTED_DUE_TO_B20_VARIATION], "Number of configurations in the scan that were rejected due to the max_B20_variation_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_L_grad_grad_B", filters[REJECTED_DUE_TO_L_GRAD_GRAD_B], "Number of configurations in the scan that were rejected due to the min_L_grad_grad_B_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_d2_volume_d_psi2", filters[REJECTED_DUE_TO_D2_VOLUME_D_PSI2], "Number of configurations in the scan that were rejected due to the max_d2_volume_d_psi2_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_DMerc", filters[REJECTED_DUE_TO_DMERC], "Number of configurations in the scan that were rejected due to the min_DMerc_to_keep filter", "dimensionless");
  nc.put("rejected_due_to_r_singularity", filters[REJECTED_DUE_TO_R_SINGULARITY], "Number of configurations in the scan that were rejected due to the min_r_singularity_to_keep filter", "dimensionless");

  int keep_all_int = (int) keep_all;
  nc.put("keep_all", keep_all_int, "1 if all configurations from the scan were saved, 0 if some configurations were filtered out", "dimensionless");
  int deterministic_int = (int) deterministic;
  nc.put("deterministic", deterministic_int, "1 if a deterministic pseudo-random number generator and seed were used, 1 if a random seed based on the time was used.", "dimensionless");
  if (!keep_all) {
    nc.put("min_R0_to_keep", min_R0_to_keep, "Configurations were kept in the scan only if the major radius of the magnetic axis was at least this value", "meter");
    nc.put("min_iota_to_keep", min_iota_to_keep, "Configurations were kept in the scan only if the absolute value of the on-axis rotational transform was at least this value", "dimensionless");
    nc.put("max_elongation_to_keep", max_elongation_to_keep, "Configurations were kept in the scan only if the elongation (in the plane perpendicular to the magnetic axis) was no greater than this value at all toroidal angles", "dimensionless");
    nc.put("min_L_grad_B_to_keep", min_L_grad_B_to_keep, "Configurations were kept in the scan only if the scale length L_grad_B (eq (3.1) in Landreman J Plasma Physics (2021)) is at least this value at each toroidal angle", "meter");
    if (q.at_least_order_r2) {
      nc.put("min_L_grad_grad_B_to_keep", min_L_grad_grad_B_to_keep, "Configurations were kept in the scan only if the scale length L_grad_grad_B (eq (3.2) in Landreman J Plasma Physics (2021)) is at least this value at each toroidal angle", "meter");
      nc.put("max_B20_variation_to_keep", max_B20_variation_to_keep, "Configurations were kept in the scan only if the range (maximum minus minimum) of B20 over toroidal angle is no more than this value", "Tesla/(meter^2)");
      nc.put("min_r_singularity_to_keep", min_r_singularity_to_keep, "Configurations were kept in the scan only if r_singularity_robust is at least this value. r_singularity_robust is the robust estimate of the minor radius at which the flux surface shapes become singular, r_c, as detailed in section 4.2 of Landreman, J Plasma Physics (2021)", "meter");
      nc.put("max_d2_volume_d_psi2_to_keep", max_d2_volume_d_psi2_to_keep, "Configurations were kept in the scan only if the magnetic well d2_volume_d_psi2 is no more than this value. d2_volume_d_psi2 is the second derivative of flux surface volume with respect to psi, where 2*pi*psi is the toroidal flux.", "Tesla^{-2} meter^{-1}");
      nc.put("min_DMerc_to_keep", min_DMerc_to_keep, "Configurations were kept in the scan only if the Mercier stability criterion DMerc_times_r2 is at least this value. DMerc_times_r2 corresponds to r^2 times the quantity DMerc in Landreman and Jorge, J Plasma Phys (2020).", "Tesla^{-2} meter^{-2}");
    }
  }
  nc.put("max_newton_iterations", q.max_newton_iterations, "Maximum iterations of Newton's method for solving the sigma equation", "dimensionless");
  nc.put("max_linesearch_iterations", q.max_linesearch_iterations, "Maximum number of times the step size is reduced in the line search for each iteration of Newton's method when solving the sigma equation", "dimensionless");
  nc.put("newton_tolerance", q.newton_tolerance, "L2 norm of the residual used as a stopping criterion for Newton's method when solving the sigma equation", "dimensionless");
  nc.put("I2", q.I2, "r^2 term in I(r), which is the toroidal current inside the flux surface times mu0/(2pi)", "Tesla/meter");
  nc.put("d_phi", q.d_phi, "Grid spacing in phi", "dimensionless");
  nc.put("B0", q.B0, "Magnetic field magnitude on the magnetic axis", "Telsa");
  nc.put("sG", q.sG, "Sign of G0", "dimensionless");
  nc.put("spsi", q.spsi, "Sign of the toroidal flux psi", "dimensionless");
  if (q.at_least_order_r2) {
    nc.put("p2", q.p2, "r^2 term in p(r), the pressure profile", "Pascal/(meter^2)");
  }
  /*
  //nc.put("axis_nmax_plus_1", R0c.size(), "Length of the arrays R0c, Z0s, etc", "dimensionless");
  nc.put("axis_length", axis_length, "Total length of the magnetic axis, from phi = 0 to 2pi", "meter");
  nc.put("d_l_d_varphi", d_l_d_varphi, "Differential arclength of the magnetic axis with respect to the Boozer toroidal angle", "meter");
  nc.put("B0_over_abs_G0", B0_over_abs_G0, "", "1/meter");
  nc.put("abs_G0_over_B0", abs_G0_over_B0, "", "meter");
  nc.put("rms_curvature", rms_curvature, "Root-mean-squared curvature of the magnetic axis, where the average is taken with respect to arclength", "1/meter");
  nc.put("grid_max_curvature", grid_max_curvature, "Maximum curvature of the magnetic axis, maximizing only over the phi grid points and not interpolating in between", "1/meter");
  nc.put("grid_max_elongation", grid_max_elongation, "Maximum elongation (ratio of major to minor axes of the O(r^1) elliptical surfaces in the plane perpendicular to the magnetic axis), maximizing only over the phi grid points and not interpolating in between", "dimensionless");
  nc.put("grid_min_R0", grid_min_R0, "Minimum major radius of the magnetic axis, minimizing only over the phi grid points and not interpolating in between", "meter");
  nc.put("grid_min_L_grad_B", grid_min_L_grad_B, "Minimum of L_grad_B over the phi grid points", "meter");
  nc.put("mean_elongation", mean_elongation, "Average elongation (ratio of major to minor axes of the O(r^1) elliptical surfaces in the plane perpendicular to the magnetic axis), where the average is taken with respect to arclength", "dimensionless");
  nc.put("mean_R", mean_R, "Average major radius of the magnetic axis, where the average is taken with respect to arclength", "meter");
  nc.put("mean_Z", mean_Z, "Average Z coordinate of the magnetic axis, where the average is taken with respect to arclength", "meter");
  nc.put("standard_deviation_of_R", standard_deviation_of_R, "Standard deviation of the major radius of the magnetic axis, where the average is taken with respect to arclength", "meter");
  nc.put("standard_deviation_of_Z", standard_deviation_of_Z, "Standard deviation of the Z coordinate of the magnetic axis, where the average is taken with respect to arclength", "meter");
  nc.put("iota", iota, "Rotational transform", "dimensionless");
  nc.put("iota_N", iota_N, "Rotational transform minus N", "dimensionless");
  if (at_least_order_r2) {
    nc.put("G2", G2, "r^2 term in G(r), which is the poloidal current outside the flux surface times mu0/(2pi)", "Tesla/meter");
    nc.put("beta_1s", beta_1s, "r * sin(theta) component of beta, the coefficient of grad psi in the Boozer covariant representation of B", "meter^{-2}");
    nc.put("B20_mean", B20_mean, "", "Tesla/(meter^2)");
    nc.put("B20_residual", B20_residual, "", "Telsa/(meter^2)");
    nc.put("B20_grid_variation", B20_grid_variation, "", "Telsa/(meter^2)");
    nc.put("d2_volume_d_psi2", d2_volume_d_psi2, "", "");
    nc.put("DGeod_times_r2", DGeod_times_r2, "", "");
    nc.put("DWell_times_r2", DWell_times_r2, "", "");
    nc.put("DMerc_times_r2", DMerc_times_r2, "", "");
    nc.put("grid_min_L_grad_grad_B", grid_min_L_grad_grad_B, "Minimum of L_grad_grad_B over the phi grid points", "meter");
    nc.put("r_singularity_robust", r_singularity_robust, "Robust estimate of the minor radius at which the flux surface shapes become singular, r_c, as detailed in section 4.2 of Landreman, J Plasma Physics (2021)", "meter");
  }
  */

  // 1D arrays
  nc.put(nphi_dim, "phi", q.phi, "The grid in the standard toroidal angle phi", "dimensionless");
  nc.put(n_scan_dim, "scan_eta_bar", scan_eta_bar, "For each configuration kept from the scan, the constant equal to B1c / B0", "1/meter");
  nc.put(n_scan_dim, "scan_sigma0", scan_sigma0, "For each configuration kept from the scan, the value of sigma at phi=0", "dimensionless");
  if (q.at_least_order_r2) {
    nc.put(n_scan_dim, "scan_B2c", scan_B2c, "For each configuration kept from the scan, the r^2 * cos(2*theta) term in |B|", "Tesla/(meter^2)");
    nc.put(n_scan_dim, "scan_B2s", scan_B2s, "For each configuration kept from the scan, the r^2 * sin(2*theta) term in |B|", "Tesla/(meter^2)");
  }
  nc.put(n_scan_dim, "scan_iota", scan_iota, "For each configuration kept from the scan, the rotational transform on axis", "dimensionless");
  nc.put(n_scan_dim, "scan_min_R0", scan_min_R0, "For each configuration kept from the scan, the minimum value of R0, the major radius of the magnetic axis. This variable corresponds to grid_min_R0 in a single Qsc calculation.", "meter");
  nc.put(n_scan_dim, "scan_max_curvature", scan_max_curvature, "For each configuration kept from the scan, the maximum curvature of the magnetic axis", "1/meter");
  nc.put(n_scan_dim, "scan_max_elongation", scan_max_elongation, "For each configuration kept from the scan, the maximum along the magnetic axis of the elongation in the plane perpendicular to the axis", "dimensionless");
  nc.put(n_scan_dim, "scan_min_L_grad_B", scan_min_L_grad_B, "For each configuration kept from the scan, the minimum along the magnetic axis of the scale length L_grad_B, (eq (3.1) in Landreman J Plasma Physics (2021). This quantity corresponds to grid_min_L_grad_B for a single Qsc run.", "meter");
  nc.put(n_scan_dim, "scan_helicity", scan_helicity, "For each configuration kept from the scan, the number of times the normal vector of the magnetic axis rotates poloidally as the axis is followed toroidally for one field period. The integer N appearing in our papers is equal to -helicity * nfp.", "dimensionless");
  nc.put(n_scan_dim, "scan_standard_deviation_of_R", scan_standard_deviation_of_R, "Standard deviation of the major radius of the magnetic axis, with respect to arclength along the axis", "meter");
  nc.put(n_scan_dim, "scan_standard_deviation_of_Z", scan_standard_deviation_of_Z, "Standard deviation of the Cartesian Z coordinate of the magnetic axis, with respect to arclength along the axis", "meter");
  if (q.at_least_order_r2) {
    nc.put(n_scan_dim, "scan_min_L_grad_grad_B", scan_min_L_grad_grad_B, "For each configuration kept from the scan, the minimum along the magnetic axis of the scale length L_grad_grad_B, (eq (3.2) in Landreman J Plasma Physics (2021). This quantity corresponds to grid_min_L_grad_grad_B for a single Qsc run.", "meter");
    nc.put(n_scan_dim, "scan_B20_variation", scan_B20_variation, "For each configuration kept from the scan, the maximum of B20 along the magnetic axis minus the minimum of B20. This quantity corresponds to B20_grid_variation for a single Qsc run.", "Telsa/(meter^2)");
    nc.put(n_scan_dim, "scan_B20_residual", scan_B20_residual, "", "");
    nc.put(n_scan_dim, "scan_r_singularity", scan_r_singularity, "For each configuration kept from the scan, the value of r_singularity_robust. r_singularity_robust is the robust estimate of the minor radius at which the flux surface shapes become singular, r_c, as detailed in section 4.2 of Landreman, J Plasma Physics (2021)", "meter");
    nc.put(n_scan_dim, "scan_d2_volume_d_psi2", scan_d2_volume_d_psi2, "For each configuration kept from the scan, the value of magnetic well d2_volume_d_psi2, the second derivative of flux surface volume with respect to psi, where 2*pi*psi is the toroidal flux.", "Tesla^{-2} meter^{-1}");
    nc.put(n_scan_dim, "scan_DMerc_times_r2", scan_DMerc_times_r2, "For each configuration kept from the scan, the overall Mercier stability criterion times the square of the effective minor radius r. This quantity corresponds to DMerc_times_r2 for a single Qsc run. DMerc (without the r^2) corresponds to the quantity DMerc in VMEC, and to DMerc in Landreman and Jorge, J Plasma Phys (2020).", "Tesla^{-2} meter^{-2}");

  }
  /*
  nc.put(nphi_dim, "curvature", curvature, "Curvature kappa of the magnetic axis", "1/meter");
  nc.put(nphi_dim, "torsion", torsion, "Torsion tau of the magnetic axis", "1/meter");
  nc.put(nphi_dim, "sigma", sigma, "Y1c / Y1s, appearing in eq (2.14) of Landreman and Sengupta, J Plasma Physics (2019)", "dimensionless");
  nc.put(nphi_dim, "X1c", X1c, "r*cos(theta) term in X, the component of the position vector in the direction of the normal vector", "dimensionless");
  nc.put(nphi_dim, "Y1s", Y1s, "r*sin(theta) term in Y, the component of the position vector in the direction of the binormal vector", "dimensionless");
  nc.put(nphi_dim, "Y1c", Y1c, "r*cos(theta) term in Y, the component of the position vector in the direction of the binormal vector", "dimensionless");
  nc.put(axis_nmax_plus_1_dim, "R0c", R0c, "Fourier cosine(n*phi) amplitudes defining the major radius R of the magnetic axis shape", "meter");
  nc.put(axis_nmax_plus_1_dim, "R0s", R0c, "Fourier sine(n*phi) amplitudes defining the major radius R of the magnetic axis shape", "meter");
  nc.put(axis_nmax_plus_1_dim, "Z0c", Z0c, "Fourier cosine(n*phi) amplitudes defining the Z coordinate of the magnetic axis shape", "meter");
  nc.put(axis_nmax_plus_1_dim, "Z0s", Z0c, "Fourier sine(n*phi) amplitudes defining the Z coordinate of the magnetic axis shape", "meter");
  nc.put(nphi_dim, "R0", R0, "Major radius of the magnetic axis", "meter");
  nc.put(nphi_dim, "Z0", Z0, "Z coordinate of the magnetic axis", "meter");
  nc.put(nphi_dim, "R0p", R0p, "d / d phi derivative of R0", "meter");
  nc.put(nphi_dim, "Z0p", Z0p, "d / d phi derivative of Z0", "meter");
  nc.put(nphi_dim, "R0pp", R0pp, "d^2 / d phi^2 derivative of R0", "meter");
  nc.put(nphi_dim, "Z0pp", Z0pp, "d^2 / d phi^2 derivative of Z0", "meter");
  nc.put(nphi_dim, "R0ppp", R0ppp, "d^3 / d phi^3 derivative of R0", "meter");
  nc.put(nphi_dim, "Z0ppp", Z0ppp, "d^3 / d phi^3 derivative of Z0", "meter");
  nc.put(nphi_dim, "d_l_d_phi", d_l_d_phi, "Differential arclength of the magnetic axis with respect to the standard toroidal angle phi", "meter");
  nc.put(nphi_dim, "d2_l_d_phi2", d2_l_d_phi2, "Second derivative of arclength of the magnetic axis with respect to the standard toroidal angle phi", "meter");
  nc.put(nphi_dim, "elongation", elongation, "Ratio of major to minor axes of the O(r^1) elliptical surfaces in the plane perpendicular to the magnetic axis", "dimensionless");
  nc.put(nphi_dim, "Boozer_toroidal_angle", Boozer_toroidal_angle, "Boozer toroidal angle varphi", "dimensionless");
  nc.put(nphi_dim, "L_grad_B", L_grad_B, "Scale length associated with first derivatives of the magnetic field, eq (3.1) in Landreman J Plasma Physics (2021)", "meter");
  nc.put(nphi_dim, "L_grad_B_inverse", L_grad_B_inverse, "1 / L_grad_B", "1/meter");
  nc.put(nphi_dim, "d_X1c_d_varphi", d_X1c_d_varphi, "Derivative of X1c with respect to the Boozer toroidal angle varphi", "dimensionless");
  nc.put(nphi_dim, "d_Y1c_d_varphi", d_Y1c_d_varphi, "Derivative of Y1c with respect to the Boozer toroidal angle varphi", "dimensionless");
  nc.put(nphi_dim, "d_Y1s_d_varphi", d_Y1s_d_varphi, "Derivative of Y1s with respect to the Boozer toroidal angle varphi", "dimensionless");
  if (at_least_order_r2) {
    nc.put(nphi_dim, "X20", X20, "r^2*cos(0*theta) term in X, the component of the position vector in the direction of the normal vector", "1/meter");
    nc.put(nphi_dim, "X2s", X2s, "r^2*sin(2*theta) term in X, the component of the position vector in the direction of the normal vector", "1/meter");
    nc.put(nphi_dim, "X2c", X2c, "r^2*cos(2*theta) term in X, the component of the position vector in the direction of the normal vector", "1/meter");
    nc.put(nphi_dim, "Y20", Y20, "r^2*cos(0*theta) term in Y, the component of the position vector in the direction of the binormal vector", "1/meter");
    nc.put(nphi_dim, "Y2s", Y2s, "r^2*sin(2*theta) term in Y, the component of the position vector in the direction of the binormal vector", "1/meter");
    nc.put(nphi_dim, "Y2c", Y2c, "r^2*cos(2*theta) term in Y, the component of the position vector in the direction of the binormal vector", "1/meter");
    nc.put(nphi_dim, "Z20", Z20, "r^2*cos(0*theta) term in Z, the component of the position vector in the direction of the tangent vector", "1/meter");
    nc.put(nphi_dim, "Z2s", Z2s, "r^2*sin(2*theta) term in Z, the component of the position vector in the direction of the tangent vector", "1/meter");
    nc.put(nphi_dim, "Z2c", Z2c, "r^2*cos(2*theta) term in Z, the component of the position vector in the direction of the tangent vector", "1/meter");

    nc.put(nphi_dim, "B20", B20, "r^2*cos(0*theta) term in the magnetic field magnitude B", "Telsa/(meter^2)");
    nc.put(nphi_dim, "B20_anomaly", B20_anomaly, "B20 - B20_mean, i.e. the toroidal variation of B that breaks O(r^2) quasisymmetry", "Telsa/(meter^2)");
	   
    nc.put(nphi_dim, "d_X20_d_varphi", d_X20_d_varphi, "Derivative of X20 with respect to the Boozer toroidal angle varphi", "1/meter");
    nc.put(nphi_dim, "d_X2s_d_varphi", d_X2s_d_varphi, "Derivative of X2s with respect to the Boozer toroidal angle varphi", "1/meter");
    nc.put(nphi_dim, "d_X2c_d_varphi", d_X2c_d_varphi, "Derivative of X2c with respect to the Boozer toroidal angle varphi", "1/meter");
    nc.put(nphi_dim, "d_Y20_d_varphi", d_Y20_d_varphi, "Derivative of Y20 with respect to the Boozer toroidal angle varphi", "1/meter");
    nc.put(nphi_dim, "d_Y2s_d_varphi", d_Y2s_d_varphi, "Derivative of Y2s with respect to the Boozer toroidal angle varphi", "1/meter");
    nc.put(nphi_dim, "d_Y2c_d_varphi", d_Y2c_d_varphi, "Derivative of Y2c with respect to the Boozer toroidal angle varphi", "1/meter");
    nc.put(nphi_dim, "d_Z20_d_varphi", d_Z20_d_varphi, "Derivative of Z20 with respect to the Boozer toroidal angle varphi", "1/meter");
    nc.put(nphi_dim, "d_Z2s_d_varphi", d_Z2s_d_varphi, "Derivative of Z2s with respect to the Boozer toroidal angle varphi", "1/meter");
    nc.put(nphi_dim, "d_Z2c_d_varphi", d_Z2c_d_varphi, "Derivative of Z2c with respect to the Boozer toroidal angle varphi", "1/meter");
    
    nc.put(nphi_dim, "L_grad_grad_B", L_grad_grad_B, "Scale length associated with second derivatives of the magnetic field, eq (3.2) in Landreman J Plasma Physics (2021)", "meter");
    nc.put(nphi_dim, "L_grad_grad_B_inverse", L_grad_grad_B_inverse, "1 / L_grad_grad_B", "1/meter");
    nc.put(nphi_dim, "r_hat_singularity_robust", r_hat_singularity_robust, "Robust estimate of the minor radius at which the flux surface shapes become singular, hat{r}_c(varphi), as detailed in section 4.2 of Landreman, J Plasma Physics (2021)", "meter");
  }
  */
  /*
  nc.put(nphi_dim, "", , "", "");
  */  
  
  // Done defining the NetCDF data.
  nc.write_and_close();
  
  if (verbose > 0) {
    end_time = std::clock();
    auto end = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for write_netcdf from chrono:           "
              << elapsed.count() << " seconds" << std::endl;
    std::cout << "Time for write_netcdf from ctime (CPU time): "
              << double(end_time - start_time) / CLOCKS_PER_SEC
              << " seconds" << std::endl;
  }
}
