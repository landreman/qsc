#include <vector>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include <netcdf.h>
#include "qsc.hpp"
#include "netcdf_writer.hpp"

using namespace qsc;

qsc::NetCDFWriter::NetCDFWriter(std::string filename) {
  int retval;
  // If you want to store long longs, "|NC_NETCDF4" must be included
  // in the file type, but such files cannot be read with
  // scipy.io.netcdf.
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

void qsc::NetCDFWriter::put(std::string varname, big& val, std::string att, std::string units) {
  // Variant for "bigs".
  // Convert results to qscfloat, since long ints require netcdf4, which scipy.io.netcdf cannot read
  qscfloat* floatval = new qscfloat;
  *floatval = (qscfloat) val;
  int var_id, retval;
  //if ((retval = nc_def_var(ncid, varname.c_str(), NC_UINT64, 0, NULL, &var_id)))
  if ((retval = nc_def_var(ncid, varname.c_str(), QSCFLOAT, 0, NULL, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(QSC_NC_FLOAT);
  pointers.push_back((void*) floatval);
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
    } else if (types[j] == QSC_NC_BIG) {
      // bigs
      std::cout << "About to nc_put a big" << std::endl;
      if ((retval = nc_put_var_ulonglong(ncid, var_ids[j], (big*) pointers[j])))
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

