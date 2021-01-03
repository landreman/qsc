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
    enum {QSC_NC_INT, QSC_NC_FLOAT, QSC_NC_BIG, QSC_NC_STRING};
    std::vector<int> types;
    static void ERR(int);
    
  public:
    NetCDFWriter(std::string);
    dim_id_type dim(std::string, int);
    void add_attribute(int, std::string, std::string);
    
    // Scalars:
    void put(std::string, int&, std::string, std::string);
    void put(std::string, qscfloat&, std::string, std::string);
    void put(std::string, big&, std::string, std::string);
    void put(std::string, std::string&, std::string);
    // 1D vectors:
    void put(dim_id_type, std::string, std::valarray<int>&, std::string, std::string);
    void put(dim_id_type, std::string, Vector&, std::string, std::string);
    // ND vectors for N > 1:
    void put(std::vector<dim_id_type>, std::string, qscfloat*, std::string, std::string);
    
    void write_and_close();
  };
}
