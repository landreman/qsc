#ifndef QSC_TOML_UTIL_H
#define QSC_TOML_UTIL_H

#include <vector>
#include <valarray>
#include "qsc.hpp"

namespace qsc {

  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, int& var);
  
  /** Handle bools
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, bool& var);

  /** Handle floats
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, qscfloat& var);

  /** Handle strings
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, std::string& var);

  /** Handle Vectors
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, Vector& var);
  
  /** Handle valarray<bool>
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, std::valarray<bool>& var);
  
  /** Handle valarray<int>
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, std::valarray<int>& var);
  
  /** Handle vector<string>
   */
  void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, std::vector<std::string>& var);
  
  /** Expand a Vector to a longer size, padding with zeros.
   */
  void pad_vector(Vector& v, std::size_t newsize);

  void toml_unused(std::vector<std::string>, toml::value);
}

#endif
