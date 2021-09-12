#include <chrono>
#include <vector>
#include "qsc.hpp"
#include "netcdf_writer.hpp"

using namespace qsc;

void Qsc::write_netcdf(std::string filename) {
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) start = std::chrono::steady_clock::now();

  if (verbose > 0) std::cout << "Writing output to " << filename << std::endl;
  qsc::NetCDFWriter nc(filename, false);

  // Define dimensions
  dim_id_type nphi_dim, axis_nmax_plus_1_dim, nbt_dim;
  nphi_dim = nc.dim("nphi", nphi);
  axis_nmax_plus_1_dim = nc.dim("axis_nmax_plus_1", R0c.size());
  nbt_dim = nc.dim("normal_binormal_tangent", 3);
  
  // Scalars
  std::string general_option = "single";
  nc.put("general_option", general_option, "Whether this job was a single configuration vs a scan");
  
  int at_least_order_r2_int = (int) at_least_order_r2;
  nc.put("at_least_order_r2", at_least_order_r2_int, "1 if the O(r^2) equations were solved, 0 if not", "dimensionless");
  int order_r2p1_int = (int) order_r2p1;
  nc.put("order_r2.1", order_r2p1_int, "1 if equations (3.12) and (3.14)-(3.15) in Landreman and Sengupta (2019) were used to compute X3c1, Y3c1, and Y3s1, 0 if not", "dimensionless");
  int order_r3_int = (int) order_r3;
  nc.put("order_r3", order_r3_int, "1 if the arrays X3 and Y3 are present in this file, 0 if not", "dimensionless");

  nc.put("order_r_option", order_r_option, "Whether the Garren-Boozer equations were solved to 1st or 2nd order in the effective minor radius r");
  nc.put("nfp", nfp, "Number of field periods", "dimensionless");
  nc.put("nphi", nphi, "Number of grid points in the toroidal angle phi", "dimensionless");
  //nc.put("axis_nmax_plus_1", R0c.size(), "Length of the arrays R0c, Z0s, etc", "dimensionless");
  nc.put("eta_bar", eta_bar, "Constant equal to B1c / B0", "1/meter");
  nc.put("sigma0", sigma0, "Value of sigma at phi=0", "dimensionless");
  nc.put("I2", I2, "r^2 term in I(r), which is the toroidal current inside the flux surface times mu0/(2pi)", "Tesla/meter");
  nc.put("d_phi", d_phi, "Grid spacing in phi", "dimensionless");
  nc.put("B0", B0, "Magnetic field magnitude on the magnetic axis", "Telsa");
  nc.put("G0", G0, "Value on the magnetic axis of G(r), which is the poloidal current outside the flux surface times mu0/(2pi)", "Tesla*meter");
  nc.put("sG", sG, "Sign of G0", "dimensionless");
  nc.put("spsi", spsi, "Sign of the toroidal flux psi", "dimensionless");
  nc.put("axis_length", axis_length, "Total length of the magnetic axis, from phi = 0 to 2pi", "meter");
  nc.put("d_l_d_varphi", d_l_d_varphi, "Differential arclength of the magnetic axis with respect to the Boozer toroidal angle", "meter");
  nc.put("B0_over_abs_G0", B0_over_abs_G0, "", "1/meter");
  nc.put("abs_G0_over_B0", abs_G0_over_B0, "", "meter");
  nc.put("helicity", helicity, "Number of times the normal vector of the magnetic axis rotates poloidally as the axis is followed toroidally for one field period. The integer N appearing in our papers is equal to -helicity * nfp.", "dimensionless");
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
  nc.put("max_newton_iterations", max_newton_iterations, "Maximum iterations of Newton's method for solving the sigma equation", "dimensionless");
  nc.put("max_linesearch_iterations", max_linesearch_iterations, "Maximum number of times the step size is reduced in the line search for each iteration of Newton's method when solving the sigma equation", "dimensionless");
  nc.put("newton_tolerance", newton_tolerance, "L2 norm of the residual used as a stopping criterion for Newton's method when solving the sigma equation", "dimensionless");
  nc.put("iota", iota, "Rotational transform", "dimensionless");
  nc.put("iota_N", iota_N, "Rotational transform minus N", "dimensionless");
  if (at_least_order_r2) {
    nc.put("B2c", B2c, "r^2 * cos(2*theta) term in |B|", "Tesla/(meter^2)");
    nc.put("B2s", B2s, "r^2 * sin(2*theta) term in |B|", "Tesla/(meter^2)");
    nc.put("p2", p2, "r^2 term in p(r), the pressure profile", "Pascal/(meter^2)");
    nc.put("G2", G2, "r^2 term in G(r), which is the poloidal current outside the flux surface times mu0/(2pi)", "Tesla/meter");
    nc.put("beta_1s", beta_1s, "r * sin(theta) component of beta, the coefficient of grad psi in the Boozer covariant representation of B", "meter^{-2}");
    nc.put("B20_mean", B20_mean, "Average over arclength along the magnetic axis of B20", "Tesla/(meter^2)");
    nc.put("B20_residual", B20_residual, "", "Telsa/(meter^2)");
    nc.put("B20_grid_variation", B20_grid_variation, "Maximum of B20 along the magnetic axis minus the minimum of B20 along the magnetic axis. The maximum and minimum are evaluated on the phi grid, without interpolation.", "Telsa/(meter^2)");
    nc.put("d2_volume_d_psi2", d2_volume_d_psi2, "Magnetic well, the second derivative of flux surface volume with respect to psi, where 2*pi*psi is the toroidal flux. Negative values are favorable for stability.", "Tesla^{-2} meter^{-1}");
    nc.put("DGeod_times_r2", DGeod_times_r2, "Geodesic curvature term in Mercier's criterion, times the square of the effective minor radius r. DGeod (without the r^2) corresponds to the quantity DGeod in VMEC, and to DGeod in Landreman and Jorge, J Plasma Phys (2020).", "Tesla^{-2} meter^{-2}");
    nc.put("DWell_times_r2", DWell_times_r2, "Magnetic well term in Mercier's criterion, times the square of the effective minor radius r. DWell (without the r^2) corresponds to the quantity DWell in VMEC, and to DWell in Landreman and Jorge, J Plasma Phys (2020).", "Tesla^{-2} meter^{-2}");
    nc.put("DMerc_times_r2", DMerc_times_r2, "Overall Mercier stability criterion times the square of the effective minor radius r. DMerc (without the r^2) corresponds to the quantity DMerc in VMEC, and to DMerc in Landreman and Jorge, J Plasma Phys (2020).", "Tesla^{-2} meter^{-2}");
    nc.put("grid_min_L_grad_grad_B", grid_min_L_grad_grad_B, "Minimum of L_grad_grad_B over the phi grid points", "meter");
    nc.put("r_singularity_robust", r_singularity_robust, "Robust estimate of the minor radius at which the flux surface shapes become singular, r_c, as detailed in section 4.2 of Landreman, J Plasma Physics (2021)", "meter");
  }
  /*
  nc.put("", , "", "");
  */

  // 1D arrays
  nc.put(nphi_dim, "phi", phi, "The grid in the standard toroidal angle phi", "dimensionless");
  nc.put(nphi_dim, "curvature", curvature, "Curvature kappa of the magnetic axis", "1/meter");
  nc.put(nphi_dim, "torsion", torsion, "Torsion tau of the magnetic axis", "1/meter");
  nc.put(nphi_dim, "sigma", sigma, "Y1c / Y1s, appearing in eq (2.14) of Landreman and Sengupta, J Plasma Physics (2019)", "dimensionless");
  nc.put(nphi_dim, "X1c", X1c, "r*cos(theta) term in X, the component of the position vector in the direction of the normal vector", "dimensionless");
  nc.put(nphi_dim, "Y1s", Y1s, "r*sin(theta) term in Y, the component of the position vector in the direction of the binormal vector", "dimensionless");
  nc.put(nphi_dim, "Y1c", Y1c, "r*cos(theta) term in Y, the component of the position vector in the direction of the binormal vector", "dimensionless");
  nc.put(axis_nmax_plus_1_dim, "R0c", R0c, "Fourier cosine(n*phi) amplitudes defining the major radius R of the magnetic axis shape", "meter");
  nc.put(axis_nmax_plus_1_dim, "R0s", R0s, "Fourier sine(n*phi) amplitudes defining the major radius R of the magnetic axis shape", "meter");
  nc.put(axis_nmax_plus_1_dim, "Z0c", Z0c, "Fourier cosine(n*phi) amplitudes defining the Z coordinate of the magnetic axis shape", "meter");
  nc.put(axis_nmax_plus_1_dim, "Z0s", Z0s, "Fourier sine(n*phi) amplitudes defining the Z coordinate of the magnetic axis shape", "meter");
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

  if (order_r2p1) {
    nc.put(nphi_dim, "lambda_for_XY3", lambda_for_XY3, "lambda in eq (3.15) and (3.12) of Landreman and Sengupta (2019), used to compute X3 and Y3", "1/meter^2");
  }

  if (order_r3) {
    nc.put(nphi_dim, "X3s1", X3s1, "r^3*sin(1*theta) term in X, the component of the position vector in the direction of the normal vector", "1/meter^2");
    nc.put(nphi_dim, "X3s3", X3s3, "r^3*sin(3*theta) term in X, the component of the position vector in the direction of the normal vector", "1/meter^2");
    nc.put(nphi_dim, "X3c1", X3c1, "r^3*cos(1*theta) term in X, the component of the position vector in the direction of the normal vector", "1/meter^2");
    nc.put(nphi_dim, "X3c3", X3c3, "r^3*cos(3*theta) term in X, the component of the position vector in the direction of the normal vector", "1/meter^2");

    nc.put(nphi_dim, "Y3s1", Y3s1, "r^3*sin(1*theta) term in Y, the component of the position vector in the direction of the binormal vector", "1/meter^2");
    nc.put(nphi_dim, "Y3s3", Y3s3, "r^3*sin(3*theta) term in Y, the component of the position vector in the direction of the binormal vector", "1/meter^2");
    nc.put(nphi_dim, "Y3c1", Y3c1, "r^3*cos(1*theta) term in Y, the component of the position vector in the direction of the binormal vector", "1/meter^2");
    nc.put(nphi_dim, "Y3c3", Y3c3, "r^3*cos(3*theta) term in Y, the component of the position vector in the direction of the binormal vector", "1/meter^2");

    nc.put(nphi_dim, "Z3s1", Z3s1, "r^3*sin(1*theta) term in Z, the component of the position vector in the direction of the tangnt vector", "1/meter^2");
    nc.put(nphi_dim, "Z3s3", Z3s3, "r^3*sin(3*theta) term in Z, the component of the position vector in the direction of the tangnt vector", "1/meter^2");
    nc.put(nphi_dim, "Z3c1", Z3c1, "r^3*cos(1*theta) term in Z, the component of the position vector in the direction of the tangnt vector", "1/meter^2");
    nc.put(nphi_dim, "Z3c3", Z3c3, "r^3*cos(3*theta) term in Z, the component of the position vector in the direction of the tangnt vector", "1/meter^2");

    nc.put(nphi_dim, "d_X3s1_d_varphi", d_X3s1_d_varphi, "Derivative of X3s1 with respect to the toroidal Boozer angle varphi", "1/meter^2");
    nc.put(nphi_dim, "d_X3s3_d_varphi", d_X3s3_d_varphi, "Derivative of X3s3 with respect to the toroidal Boozer angle varphi", "1/meter^2");
    nc.put(nphi_dim, "d_X3c1_d_varphi", d_X3c1_d_varphi, "Derivative of X3c1 with respect to the toroidal Boozer angle varphi", "1/meter^2");
    nc.put(nphi_dim, "d_X3c3_d_varphi", d_X3c3_d_varphi, "Derivative of X3c3 with respect to the toroidal Boozer angle varphi", "1/meter^2");

    nc.put(nphi_dim, "d_Y3s1_d_varphi", d_Y3s1_d_varphi, "Derivative of Y3s1 with respect to the toroidal Boozer angle varphi", "1/meter^2");
    nc.put(nphi_dim, "d_Y3s3_d_varphi", d_Y3s3_d_varphi, "Derivative of Y3s3 with respect to the toroidal Boozer angle varphi", "1/meter^2");
    nc.put(nphi_dim, "d_Y3c1_d_varphi", d_Y3c1_d_varphi, "Derivative of Y3c1 with respect to the toroidal Boozer angle varphi", "1/meter^2");
    nc.put(nphi_dim, "d_Y3c3_d_varphi", d_Y3c3_d_varphi, "Derivative of Y3c3 with respect to the toroidal Boozer angle varphi", "1/meter^2");

  }
  /*
  nc.put(nphi_dim, "", , "", "");
  */
  
  // ND arrays for N > 1:
  std::vector<dim_id_type> nphi_nphi_dim {nphi_dim, nphi_dim};
  nc.put(nphi_nphi_dim, "d_d_phi", &d_d_phi(0, 0),
	 "Pseudospectral differentiation matrix with respect to the standard toroidal angle phi", "dimensionless");
  nc.put(nphi_nphi_dim, "d_d_varphi", &d_d_varphi(0, 0),
	 "Pseudospectral differentiation matrix with respect to the Boozer toroidal angle varphi", "dimensionless");

  std::vector<dim_id_type> nphi_nbt_nbt_dim {nphi_dim, nbt_dim, nbt_dim};
  nc.put(nphi_nbt_nbt_dim, "grad_B_tensor", &grad_B_tensor(0, 0, 0),
	 "The grad B tensor at each grid point along the magnetic axis, eq (3.12) in Landreman J Plasma Physics (2021)", "Tesla/meter");

  std::vector<dim_id_type> nphi_nbt_nbt_nbt_dim {nphi_dim, nbt_dim, nbt_dim, nbt_dim};
  if (at_least_order_r2) {
    nc.put(nphi_nbt_nbt_nbt_dim, "grad_grad_B_tensor", &grad_grad_B_tensor(0, 0, 0, 0),
	   "The grad grad B tensor at each grid point along the magnetic axis, eq (3.13) in Landreman J Plasma Physics (2021)", "Tesla/(meter^2)");
  }
  
  // Done defining the NetCDF data.
  nc.write_and_close();
  
  if (verbose > 0) {
    auto end = std::chrono::steady_clock::now();    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for write_netcdf: "
              << elapsed.count() << " seconds" << std::endl;
  }
}
