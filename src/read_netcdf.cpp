#include <vector>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include <netcdf.h>
#include "qsc.hpp"

using namespace qsc;

#ifdef SINGLE
#define nc_get_var_qscfloat nc_get_var_float
#else
#define nc_get_var_qscfloat nc_get_var_double
#endif

namespace qsc {
  
  /** A class to streamline the process of reading a NetCDF file.
   */
  class NetCDFReader {
  private:
    int ncid, ndims, nvars, ngatts, unlimdimid;
    static void ERR(int);
    
  public:
    NetCDFReader(std::string);
    
    // Scalars:
    void get(std::string, int&);
    void get(std::string, qscfloat&);
    // Vectors
    void get(std::string, Vector&);
    // Strings
    void get(std::string, std::string&);
	     
    void close();
  };
}

qsc::NetCDFReader::NetCDFReader(std::string filename) {
  int retval;
  if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
    ERR(retval);
  // Get the number of dimensions, variables, etc:
  if ((retval = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
      ERR(retval);
}

void qsc::NetCDFReader::ERR(int e) {
  std::cout << "NetCDF Error! " << nc_strerror(e) << std::endl;
  throw std::runtime_error(nc_strerror(e));
}

void qsc::NetCDFReader::get(std::string varname, int& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_int(ncid, var_id, &var)))
      ERR(retval);
}

void qsc::NetCDFReader::get(std::string varname, qscfloat& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_qscfloat(ncid, var_id, &var)))
      ERR(retval);
}

void qsc::NetCDFReader::get(std::string varname, Vector& var) {
  int var_id, retval;
  std::cout << "Reading field " << varname << std::endl;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_qscfloat(ncid, var_id, &var[0])))
      ERR(retval);
}

void qsc::NetCDFReader::get(std::string varname, std::string& var) {
  int var_id, retval;
  char nc_string[50]; // 50 is hard-coded as the string length in quasisymmetry_variables.f90
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_text(ncid, var_id, nc_string)))
      ERR(retval);
  // Convert char array to std::string
  //std::cout << "nc_string:" << nc_string << std::endl;
  var = std::string(nc_string);
  //std::cout << "var:" << var << " size=" << var.size() << " var[2]=" << ((int)var[2])
  //	    << " var[3]=" << ((int)var[3]) << std::endl;
  // NetCDF pads the strings with spaces (ascii 32) no ascii 0.
  int index = var.find(" ");
  //std::cout << "index=" << index << std::endl;
  var.resize(index);
  //std::cout << "After resize, var:" << var << " size=" << var.size() << std::endl;
}

void qsc::NetCDFReader::close() {
  int retval;
  if ((retval = nc_close(ncid))) ERR(retval);
}

//////////////////////////////////////////////////////////////////////
// End of definition of the qsc::NetCDFReader class.
//////////////////////////////////////////////////////////////////////

/** Read in quantities for a Qsc object from a NetCDF file. This
    function is used for testing.
 */
void Qsc::read_netcdf(std::string filename, char C_or_F) {
  if (C_or_F != 'C' && C_or_F != 'F') throw std::runtime_error("C_or_F must be C or F");
  int fortran = (C_or_F == 'F');
  
  if (verbose > 0) std::cout << "About to try reading netcdf file " << filename << std::endl;
  qsc::NetCDFReader nc(filename);

  if (fortran) {
    nc.get("N_phi", nphi);
    nc.get("order_r_option", order_r_option);
    std::cout << "order_r_option:" << order_r_option << std::endl;
    at_least_order_r2 = !(order_r_option.compare("r1") == 0);
    std::cout << "at_least_order_r2:" << at_least_order_r2 << std::endl;
    allocate();

    // Scalars
    nc.get("nfp", nfp);
    nc.get("sign_G", sG);
    nc.get("sign_psi", spsi);
    nc.get("eta_bar", eta_bar);
    nc.get("sigma_initial", sigma0);
    nc.get("iota", iota);
    nc.get("mean_elongation", mean_elongation);
    nc.get("rms_curvature", rms_curvature);
    nc.get("axis_length", axis_length);
    nc.get("abs_G0_over_B0", abs_G0_over_B0);
    nc.get("standard_deviation_of_R", standard_deviation_of_R);
    nc.get("standard_deviation_of_Z", standard_deviation_of_Z);
    nc.get("axis_helicity", helicity);
    nc.get("B0", B0);
    if (at_least_order_r2) {
      nc.get("p2", p2);
      nc.get("B2s", B2s);
      nc.get("B2c", B2c);
    }
    
    // Vectors
    nc.get("phi", phi);
    nc.get("Boozer_toroidal_angle", Boozer_toroidal_angle);
    nc.get("R0", R0);
    nc.get("z0", Z0);
    nc.get("curvature", curvature);
    nc.get("torsion", torsion);
    nc.get("sigma", sigma);
    nc.get("X1c", X1c);
    nc.get("Y1c", Y1c);
    nc.get("Y1s", Y1s);
    nc.get("elongation", elongation);    
    nc.get("d_l_d_phi", d_l_d_phi);
    nc.get("modBinv_sqrt_half_grad_B_colon_grad_B", L_grad_B_inverse);
    if (at_least_order_r2) {
      nc.get("X20", X20);
      nc.get("X2s", X2s);
      nc.get("X2c", X2c);
      nc.get("Y20", Y20);
      nc.get("Y2s", Y2s);
      nc.get("Y2c", Y2c);
      nc.get("Z20", Z20);
      nc.get("Z2s", Z2s);
      nc.get("Z2c", Z2c);
      nc.get("B20", B20);
    }
    
  } else {
    // Data saved by the C++ version
    
    nc.get("nphi", nphi);
    allocate();
    
    // Scalars
    int tempint;
    nc.get("at_least_order_r2", tempint);
    at_least_order_r2 = (bool) tempint;
    nc.get("nfp", nfp);
    nc.get("eta_bar", eta_bar);
    nc.get("sigma0", sigma0);
    nc.get("I2", I2);
    nc.get("d_phi", d_phi);
    nc.get("B0", B0);
    nc.get("G0", G0);
    nc.get("sG", sG);
    nc.get("spsi", spsi);
    nc.get("axis_length", axis_length);
    nc.get("d_l_d_varphi", d_l_d_varphi);
    nc.get("B0_over_abs_G0", B0_over_abs_G0);
    nc.get("abs_G0_over_B0", abs_G0_over_B0);
    nc.get("helicity", helicity);
    nc.get("rms_curvature", rms_curvature);
    nc.get("grid_max_curvature", grid_max_curvature);
    nc.get("grid_max_elongation", grid_max_elongation);
    nc.get("grid_min_R0", grid_min_R0);
    nc.get("mean_elongation", mean_elongation);
    nc.get("mean_R", mean_R);
    nc.get("mean_Z", mean_Z);
    nc.get("standard_deviation_of_R", standard_deviation_of_R);
    nc.get("standard_deviation_of_Z", standard_deviation_of_Z);
    nc.get("max_newton_iterations", max_newton_iterations);
    nc.get("max_linesearch_iterations", max_linesearch_iterations);
    nc.get("newton_tolerance", newton_tolerance);
    nc.get("iota", iota);
    nc.get("iota_N", iota_N);
    if (at_least_order_r2) {
      nc.get("B2c", B2c);
      nc.get("B2s", B2s);
      nc.get("p2", p2);
      nc.get("G2", G2);
      nc.get("beta_1s", beta_1s);
    }

    // Vectors
    nc.get("phi", phi);
    nc.get("curvature", curvature);
    nc.get("torsion", torsion);
    nc.get("sigma", sigma);
    nc.get("X1c", X1c);
    nc.get("Y1c", Y1c);
    nc.get("Y1s", Y1s);
    nc.get("R0", R0);
    nc.get("Z0", Z0);
    nc.get("d_l_d_phi", d_l_d_phi);
    nc.get("d2_l_d_phi2", d2_l_d_phi2);
    nc.get("elongation", elongation);
    nc.get("Boozer_toroidal_angle", Boozer_toroidal_angle);
    nc.get("L_grad_B", L_grad_B);
    nc.get("L_grad_B_inverse", L_grad_B_inverse);
    nc.get("d_X1c_d_varphi", d_X1c_d_varphi);
    nc.get("d_Y1c_d_varphi", d_Y1c_d_varphi);
    nc.get("d_Y1s_d_varphi", d_Y1s_d_varphi);
    if (at_least_order_r2) {
      nc.get("X20", X20);
      nc.get("X2s", X2s);
      nc.get("X2c", X2c);
      nc.get("Y20", Y20);
      nc.get("Y2s", Y2s);
      nc.get("Y2c", Y2c);
      nc.get("Z20", Z20);
      nc.get("Z2s", Z2s);
      nc.get("Z2c", Z2c);
      nc.get("B20", B20);
      nc.get("d_X20_d_varphi", d_X20_d_varphi);
      nc.get("d_X2s_d_varphi", d_X2s_d_varphi);
      nc.get("d_X2c_d_varphi", d_X2c_d_varphi);
      nc.get("d_Y20_d_varphi", d_Y20_d_varphi);
      nc.get("d_Y2s_d_varphi", d_Y2s_d_varphi);
      nc.get("d_Y2c_d_varphi", d_Y2c_d_varphi);
      nc.get("d_Z20_d_varphi", d_Z20_d_varphi);
      nc.get("d_Z2s_d_varphi", d_Z2s_d_varphi);
      nc.get("d_Z2c_d_varphi", d_Z2c_d_varphi);
    }
  }
  // nc.get("", );

  nc.close();
  
}
