#include <ctime>
#include <chrono>
//#include <netcdfcpp.h>
#include <netcdf.h>
//#include <netcdf>
#include "qsc.hpp"

using namespace qsc;
//using namespace netCDF;
//using namespace netCDF::exceptions;

void Qsc::write_output() {
  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }
  
  //NcFile dataFile("sfc_pres_temp.nc", NcFile::Replace);
  //NcFile dataFile("simple_xy.nc", NcFile::replace);
  
  int ncid, retval;
  retval = nc_create("qsc_out.foo.nc", NC_CLOBBER, &ncid);
  retval = nc_close(ncid);
  
  if (verbose > 0) {
    end_time = std::clock();
    auto end = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for write_output from chrono:           "
              << elapsed.count() << " seconds" << std::endl;
    std::cout << "Time for write_output from ctime (CPU time): "
              << double(end_time - start_time) / CLOCKS_PER_SEC
              << " seconds" << std::endl;
  }
}
