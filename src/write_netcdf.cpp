#include <ctime>
#include <chrono>
#include <vector>
#include <valarray>
#include <stdexcept>
#include <netcdf.h>
#include "qsc.hpp"

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
  
  class NetCDF {
  private:
    int ncid;
    std::vector<int> var_ids;
    std::vector<void*> pointers;
    enum {QSC_NC_INT, QSC_NC_FLOAT};
    std::vector<int> types;
    static void ERR(int);
  public:
    NetCDF(std::string);
    dim_id_type dim(std::string, int);
    void add_attribute(int, std::string);
    // Scalars:
    void put(std::string, int&, std::string);
    void put(std::string, qscfloat&, std::string);
    // 1D vectors:
    void put(std::string, dim_id_type, std::valarray<int>&, std::string);
    void put(std::string, dim_id_type, Vector&, std::string);
    // ND vectors for N > 1:
    //void put(std::string, int*, int, std::string);
    //void put(std::string, int*, qscfloat, std::string);
    void write_and_close();
  };
}

qsc::NetCDF::NetCDF(std::string filename) {
  int retval;
  if ((retval = nc_create(filename.c_str(), NC_CLOBBER, &ncid))) ERR(retval);
}

void qsc::NetCDF::ERR(int e) {
  throw std::runtime_error(nc_strerror(e));
}

int qsc::NetCDF::dim(std::string dimname, int val) {
  int dim_id, retval;
  if ((retval = nc_def_dim(ncid, dimname.c_str(), val, &dim_id)))
    ERR(retval);
  return dim_id;
}

/**
 * If an empty string is provided, no attribute is written.
 */
void qsc::NetCDF::add_attribute(int var_id, std::string str) {
  int retval;
  if (str.size() > 0) {
    if ((retval = nc_put_att_text(ncid, var_id, "description", 
				  str.size(), str.c_str())))
      ERR(retval);
  }
}

void qsc::NetCDF::put(std::string varname, int& val, std::string att) {
  // Variant for scalar ints
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), NC_INT, 0, NULL, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_INT);
  pointers.push_back((void*) &val);
  add_attribute(var_id, att);
}

void qsc::NetCDF::put(std::string varname, qscfloat& val, std::string att) {
  // Variant for scalar floats
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), QSCFLOAT, 0, NULL, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_FLOAT);
  pointers.push_back((void*) &val);
  add_attribute(var_id, att);
}

void qsc::NetCDF::put(std::string varname, dim_id_type dim_id, std::valarray<int>& val, std::string att) {
  // Variant for 1D int arrays
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), NC_INT, 1, &dim_id, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_INT);
  pointers.push_back((void*) &val[0]);
  add_attribute(var_id, att);
}

void qsc::NetCDF::put(std::string varname, dim_id_type dim_id, Vector& val, std::string att) {
  // Variant for 1D float arrays
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), QSCFLOAT, 1, &dim_id, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_FLOAT);
  pointers.push_back((void*) &val[0]);
  add_attribute(var_id, att);
}

void qsc::NetCDF::write_and_close() {
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
// End of definition of the qsc::NetCDF class.
//////////////////////////////////////////////////////////////////////

void Qsc::write_netcdf() {
  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }
  
  qsc::NetCDF nc("qsc_out.foo.nc");

  // Define dimensions
  dim_id_type nphi_dim, axis_nmax_plus_1_dim;
  nphi_dim = nc.dim("nphi", nphi);
  axis_nmax_plus_1_dim = nc.dim("axis_nmax_plus_1", R0c.size() + 1);
  
  // Scalars
  nc.put("nfp", nfp, "Number of field periods");
  nc.put("eta_bar", eta_bar, "Constant equal to B1c / B0");

  // 1D arrays
  nc.put("curvature", nphi_dim, curvature, "Curvature kappa of the magnetic axis");
  nc.put("torsion", nphi_dim, torsion, "Torsion tau of the magnetic axis");
  nc.put("sigma", nphi_dim, sigma, "");
  nc.put("X1c", nphi_dim, X1c, "");
  nc.put("Y1s", nphi_dim, Y1s, "");
  nc.put("Y1c", nphi_dim, Y1c, "");

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
