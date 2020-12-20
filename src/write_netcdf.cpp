#include <ctime>
#include <chrono>
//#include <netcdfcpp.h>
#include <netcdf.h>
//#include <netcdf>
#include <stdexcept>
#include "qsc.hpp"

using namespace qsc;
//using namespace netCDF;
//using namespace netCDF::exceptions;

#ifdef SINGLE
#define QSCFLOAT NC_FLOAT
#define nc_put_var_qscfloat nc_put_var_float
#else
#define QSCFLOAT NC_DOUBLE
#define nc_put_var_qscfloat nc_put_var_double
#endif

#define ERR(e) {throw std::runtime_error(nc_strerror(e));}

void Qsc::write_netcdf() {
  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }
  
  //NcFile dataFile("sfc_pres_temp.nc", NcFile::Replace);
  //NcFile dataFile("simple_xy.nc", NcFile::replace);

  
  int ncid, retval;
  if ((retval = nc_create("qsc_out.foo.nc", NC_CLOBBER, &ncid))) ERR(retval);

  // Define dimensions
  int nphi_dim, axis_nmax_plus_1_dim;
  if ((retval = nc_def_dim(ncid, "nphi", nphi, &nphi_dim))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "axis_nmax_plus_1", R0c.size() + 1, &axis_nmax_plus_1_dim))) ERR(retval);

  // Scalars
  int nfp_id; if ((retval = nc_def_var(ncid, "nfp", NC_INT, 0, NULL, &nfp_id))) ERR(retval);
  int eta_bar_id; if ((retval = nc_def_var(ncid, "eta_bar", QSCFLOAT, 0, NULL, &eta_bar_id))) ERR(retval);

  // 1D arrays
  int sigma_id; if ((retval = nc_def_var(ncid, "sigma", QSCFLOAT, 1, &nphi_dim, &sigma_id))) ERR(retval);
  
  // End define mode. This tells netCDF we are done defining metadata.
  if ((retval = nc_enddef(ncid))) ERR(retval);

  ///////////////////////////////////////////////////////////////////
  // Write the data
  ///////////////////////////////////////////////////////////////////

  // Scalars
  if ((retval = nc_put_var_int(ncid, nfp_id, &nfp))) ERR(retval);
  if ((retval = nc_put_var_qscfloat(ncid, eta_bar_id, &eta_bar))) ERR(retval);

  // 1D arrays
  if ((retval = nc_put_var_qscfloat(ncid, sigma_id, &sigma[0]))) ERR(retval);
  
  // Close the file
  if ((retval = nc_close(ncid))) ERR(retval);
  
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
