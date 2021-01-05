#include <iostream>
#include <stdexcept>
#include "toml.hpp"
#include "toml_util.hpp"

using namespace qsc;

void qsc::toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, int& var) {
  // toml11's find_or does not seem to give informative error messages
  // like its "find" method does, hence my little function here.
  if (indata.contains(varname)) {
    var = toml::find<int>(indata, varname);
  }
  varlist.push_back(varname);
  /*
  try {
    var = toml::find_or(indata, varname, var);
  } catch (std::exception) {
    std::cerr << "Unable to read variable " << varname << " from input file " << filename << std::endl;
    throw;
  }
  */
}

/** Handle bools
 */
void qsc::toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, bool& var) {
  if (indata.contains(varname)) {
    var = toml::find<bool>(indata, varname);
  }
  varlist.push_back(varname);
}

/** Handle floats
 */
void qsc::toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, qscfloat& var) {
  if (indata.contains(varname)) {
    // toml11 documentation suggests the following to handle the case
    // in which an int is provided rather than a float.
    auto vx = toml::find(indata, varname);
    if (vx.is_floating()) {
      var = vx.as_floating(std::nothrow);
    } else {
      var = static_cast<double>(vx.as_integer());
    }
    //var = toml::find<double>(indata, varname);
  }
  varlist.push_back(varname);
}

/** Handle strings
 */
void qsc::toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, std::string& var) {
  if (indata.contains(varname)) {
    var = toml::find<std::string>(indata, varname);
  }
  varlist.push_back(varname);
}

/** Handle Vectors
 */
void qsc::toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, Vector& var) {
  if (indata.contains(varname)) {
    auto indata_vector = toml::find<std::vector<double>>(indata, varname);
    // Convert the std::vector<double> to a Vector:
    var.resize(indata_vector.size(), 0.0);
    for (int j = 0; j < indata_vector.size(); j++) var[j] = (qscfloat)indata_vector[j];
  }
  varlist.push_back(varname);
}

/** Expand a Vector to a longer size, padding with zeros.
 */
void qsc::pad_vector(Vector& v, std::size_t newsize) {
  assert(v.size() <= newsize);
  Vector vcopy(v);

  v.resize(newsize, 0.0);
  for (int j = 0; j < vcopy.size(); j++) v[j] = vcopy[j];
}

/** Check to see if there are any unused keys in the input file.
 */
void qsc::toml_unused(std::vector<std::string> varlist, toml::value indata) {
  int j;
  /*
  std::cout << "varlist:";
  for (j = 0; j < varlist.size(); j++) std::cout << " " << varlist[j];
  std::cout << std::endl;
  */
  
  int found_match;
  for (auto item : indata.as_table()) {
    // I believe "item" has type std::pair, representing a key-value pair.
    auto key = item.first;
    //std::cout << "Found key: " << key << std::endl;
    found_match = 0;
    for (j = 0; j < varlist.size(); j++) {
      if (key.compare(varlist[j]) == 0) {
	found_match = 1;
	break;
      }
    }
    if (found_match == 0) {
      throw std::runtime_error(std::string("Unused key in input file: ") + key);
    }
  }
}
