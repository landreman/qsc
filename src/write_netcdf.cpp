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
    void add_attribute(int, std::string);
    
    // Scalars:
    void put(std::string, int&, std::string);
    void put(std::string, qscfloat&, std::string);
    // 1D vectors:
    void put(dim_id_type, std::string, std::valarray<int>&, std::string);
    void put(dim_id_type, std::string, Vector&, std::string);
    // ND vectors for N > 1:
    void put(std::vector<dim_id_type>, std::string, qscfloat*, std::string);
    
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
void qsc::NetCDFWriter::add_attribute(int var_id, std::string str) {
  int retval;
  if (str.size() > 0) {
    if ((retval = nc_put_att_text(ncid, var_id, "description", 
				  str.size(), str.c_str())))
      ERR(retval);
  }
}

void qsc::NetCDFWriter::put(std::string varname, int& val, std::string att) {
  // Variant for scalar ints
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), NC_INT, 0, NULL, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_INT);
  pointers.push_back((void*) &val);
  add_attribute(var_id, att);
}

void qsc::NetCDFWriter::put(std::string varname, qscfloat& val, std::string att) {
  // Variant for scalar floats
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), QSCFLOAT, 0, NULL, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_FLOAT);
  pointers.push_back((void*) &val);
  add_attribute(var_id, att);
}

void qsc::NetCDFWriter::put(dim_id_type dim_id, std::string varname, std::valarray<int>& val, std::string att) {
  // Variant for 1D int arrays
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), NC_INT, 1, &dim_id, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_INT);
  pointers.push_back((void*) &val[0]);
  add_attribute(var_id, att);
}

void qsc::NetCDFWriter::put(dim_id_type dim_id, std::string varname, Vector& val, std::string att) {
  // Variant for 1D float arrays
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), QSCFLOAT, 1, &dim_id, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_FLOAT);
  pointers.push_back((void*) &val[0]);
  add_attribute(var_id, att);
}

void qsc::NetCDFWriter::put(std::vector<dim_id_type> dim_id, std::string varname, qscfloat* pointer, std::string att) {
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
  add_attribute(var_id, att);
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
  axis_nmax_plus_1_dim = nc.dim("axis_nmax_plus_1", R0c.size() + 1);
  nbt_dim = nc.dim("normal_binormal_tangent", 3);
  
  // Scalars
  nc.put("nfp", nfp, "Number of field periods");
  nc.put("eta_bar", eta_bar, "Constant equal to B1c / B0");

  // 1D arrays
  nc.put(nphi_dim, "phi", phi, "The grid in the standard toroidal angle phi");
  nc.put(nphi_dim, "curvature", curvature, "Curvature kappa of the magnetic axis");
  nc.put(nphi_dim, "torsion", torsion, "Torsion tau of the magnetic axis");
  nc.put(nphi_dim, "sigma", sigma, "");
  nc.put(nphi_dim, "X1c", X1c, "");
  nc.put(nphi_dim, "Y1s", Y1s, "");
  nc.put(nphi_dim, "Y1c", Y1c, "");

  // ND arrays for N > 1:
  std::vector<dim_id_type> nphi_nphi_dim {nphi_dim, nphi_dim};
  nc.put(nphi_nphi_dim, "d_d_phi", &d_d_phi(0, 0),
	 "Pseudospectral differentiation matrix with respect to the standard toroidal angle phi");
  nc.put(nphi_nphi_dim, "d_d_varphi", &d_d_varphi(0, 0),
	 "Pseudospectral differentiation matrix with respect to the Boozer toroidal angle varphi");

  std::vector<dim_id_type> nphi_nbt_nbt_dim {nphi_dim, nbt_dim, nbt_dim};
  nc.put(nphi_nbt_nbt_dim, "grad_B_tensor", &grad_B_tensor(0, 0, 0),
	 "The grad B tensor at each grid point along the magnetic axis, eq (3.12) in Landreman J Plasma Physics (2021)");

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
