#ifndef QSC_TOML_UTIL_H
#define QSC_TOML_UTIL_H

#include <vector>
#include "qsc.hpp"

namespace qsc {

  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, int& var);
  
  /** Handle floats
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, qscfloat& var);

  /** Handle strings
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, std::string& var);

  /** Handle Vectors
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, Vector& var);
  
  /** Expand a Vector to a longer size, padding with zeros.
   */
  void pad_vector(Vector& v, std::size_t newsize);

  void toml_unused(std::vector<std::string>, toml::value);
}

#endif
