#include <ctime>
#include <chrono>
#include <vector>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include <netcdf.h>
#include "qsc.hpp"

/** QSC uses the C interface to NetCDF, not the C++ interface. One
    reason for this is to simplify the build system and avoid a
    potential problem for users, since the C++ interface does not
    automatically come with NetCDF. Second, I never managed to get the
    build system to work with the C++ interface on my laptop. Third,
    NetCDF's C++ interface is rather clunky, requiring multiple lines
    of code to add each variable, so I wanted to write a wrapper class
    to simplify this anyway.
 */

using namespace qsc;

#ifdef SINGLE
#define QSCFLOAT NC_FLOAT
#define nc_put_var_qscfloat nc_put_var_float
#else
#define QSCFLOAT NC_DOUBLE
#define nc_put_var_qscfloat nc_put_var_double
#endif

namespace qsc {
  
  typedef int dim_id_type;

  /** A class to streamline the process of writing a NetCDF file.
   */
  class NetCDFWriter {
  private:
    int ncid;
    std::vector<int> var_ids;
    std::vector<void*> pointers;
    enum {QSC_NC_INT, QSC_NC_FLOAT};
    std::vector<int> types;
    static void ERR(int);
    
  public:
    NetCDFWriter(std::string);
    dim_id_type dim(std::string, int);
    void add_attribute(int, std::string, std::string);
    
    // Scalars:
    void put(std::string, int&, std::string, std::string);
    void put(std::string, qscfloat&, std::string, std::string);
    // 1D vectors:
    void put(dim_id_type, std::string, std::valarray<int>&, std::string, std::string);
    void put(dim_id_type, std::string, Vector&, std::string, std::string);
    // ND vectors for N > 1:
    void put(std::vector<dim_id_type>, std::string, qscfloat*, std::string, std::string);
    
    void write_and_close();
  };
}

qsc::NetCDFWriter::NetCDFWriter(std::string filename) {
  int retval;
  if ((retval = nc_create(filename.c_str(), NC_CLOBBER, &ncid))) ERR(retval);
}

void qsc::NetCDFWriter::ERR(int e) {
  throw std::runtime_error(nc_strerror(e));
}

int qsc::NetCDFWriter::dim(std::string dimname, int val) {
  int dim_id, retval;
  if ((retval = nc_def_dim(ncid, dimname.c_str(), val, &dim_id)))
    ERR(retval);
  return dim_id;
}

/**
 * If an empty string is provided, no attribute is written.
 */
void qsc::NetCDFWriter::add_attribute(int var_id, std::string str, std::string units) {
  int retval;
  if (str.size() > 0) {
    if ((retval = nc_put_att_text(ncid, var_id, "description",
				  str.size(), str.c_str())))
      ERR(retval);
  }
  if (units.size() > 0) {
    if ((retval = nc_put_att_text(ncid, var_id, "units",
				  units.size(), units.c_str())))
      ERR(retval);
  }
}

void qsc::NetCDFWriter::put(std::string varname, int& val, std::string att, std::string units) {
  // Variant for scalar ints
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), NC_INT, 0, NULL, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_INT);
  pointers.push_back((void*) &val);
  add_attribute(var_id, att, units);
}

void qsc::NetCDFWriter::put(std::string varname, qscfloat& val, std::string att, std::string units) {
  // Variant for scalar floats
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), QSCFLOAT, 0, NULL, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_FLOAT);
  pointers.push_back((void*) &val);
  add_attribute(var_id, att, units);
}

void qsc::NetCDFWriter::put(dim_id_type dim_id, std::string varname, std::valarray<int>& val, std::string att, std::string units) {
  // Variant for 1D int arrays
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), NC_INT, 1, &dim_id, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_INT);
  pointers.push_back((void*) &val[0]);
  add_attribute(var_id, att, units);
}

void qsc::NetCDFWriter::put(dim_id_type dim_id, std::string varname, Vector& val, std::string att, std::string units) {
  // Variant for 1D float arrays
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), QSCFLOAT, 1, &dim_id, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_FLOAT);
  pointers.push_back((void*) &val[0]);
  add_attribute(var_id, att, units);
}

void qsc::NetCDFWriter::put(std::vector<dim_id_type> dim_id, std::string varname, qscfloat* pointer, std::string att, std::string units) {
  // Variant for ND float arrays for N > 1
  
  // NetCDF wants the order of the dimensions to be reversed compared to the QSC definitions.
  // Therefore we make a copy of the array of dimensions, and reverse the order of the copy.
  std::vector<dim_id_type> dim_id_reversed(dim_id);
  std::reverse(std::begin(dim_id_reversed), std::end(dim_id_reversed));

  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), QSCFLOAT, dim_id.size(), &dim_id_reversed[0], &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_FLOAT);
  pointers.push_back((void*) pointer);
  add_attribute(var_id, att, units);
}

void qsc::NetCDFWriter::write_and_close() {
  int retval;

  // End define mode. This tells netCDF we are done defining metadata.
  if ((retval = nc_enddef(ncid))) ERR(retval);

  // Write the data
  for (int j = 0; j < var_ids.size(); j++) {
    if (types[j] == QSC_NC_INT) {
      // ints
      if ((retval = nc_put_var_int(ncid, var_ids[j], (int*) pointers[j])))
	ERR(retval);
    } else {
      // floats
      if ((retval = nc_put_var_qscfloat(ncid, var_ids[j], (qscfloat*) pointers[j])))
	ERR(retval);
    }
  }
  
  // Close the file
  if ((retval = nc_close(ncid))) ERR(retval);
}

//////////////////////////////////////////////////////////////////////
// End of definition of the qsc::NetCDFWriter class.
//////////////////////////////////////////////////////////////////////

void Qsc::write_netcdf(std::string filename) {
  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }

  if (verbose > 0) std::cout << "Writing output to " << filename << std::endl;
  qsc::NetCDFWriter nc(filename);

  // Define dimensions
  dim_id_type nphi_dim, axis_nmax_plus_1_dim, nbt_dim;
  nphi_dim = nc.dim("nphi", nphi);
  axis_nmax_plus_1_dim = nc.dim("axis_nmax_plus_1", R0c.size());
  nbt_dim = nc.dim("normal_binormal_tangent", 3);
  
  // Scalars
  int at_least_order_r2_int = (int) at_least_order_r2;
  nc.put("at_least_order_r2", at_least_order_r2_int, "1 if the O(r^2) equations were solved, 0 if not", "dimensionless");
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
    nc.put("B20_mean", B20_mean, "", "Tesla/(meter^2)");
    nc.put("B20_residual", B20_residual, "", "Telsa/(meter^2)");
    nc.put("B20_grid_variation", B20_grid_variation, "", "Telsa/(meter^2)");
    nc.put("d2_volume_d_psi2", d2_volume_d_psi2, "", "");
    nc.put("DGeod_times_r2", DGeod_times_r2, "", "");
    nc.put("DWell_times_r2", DWell_times_r2, "", "");
    nc.put("DMerc_times_r2", DMerc_times_r2, "", "");
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
