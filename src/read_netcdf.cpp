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
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_qscfloat(ncid, var_id, &var[0])))
      ERR(retval);
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
  } else {
    nc.get("nphi", nphi);
  }
  allocate();

  // Scalars
  nc.get("nfp", nfp);
  nc.get("iota", iota);
  nc.get("mean_elongation", mean_elongation);
  nc.get("rms_curvature", rms_curvature);
  if (fortran) {
    nc.get("axis_helicity", helicity);
  } else {
    nc.get("helicity", helicity);
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
  nc.get("d_X1c_d_varphi", d_X1c_d_varphi);
  nc.get("d_Y1c_d_varphi", d_Y1c_d_varphi);
  nc.get("d_Y1s_d_varphi", d_Y1s_d_varphi);
  if (fortran) {
    nc.get("modBinv_sqrt_half_grad_B_colon_grad_B", L_grad_B_inverse);
    L_grad_B = 1.0 / L_grad_B_inverse;
  } else {
    nc.get("L_grad_B", L_grad_B);
    nc.get("L_grad_B_inverse", L_grad_B_inverse);
  }

  nc.close();
  
}
