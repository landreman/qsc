#include <exception>
#include <stdexcept>
#include <algorithm>
#include "qsc.hpp"
#include "toml.hpp"

using namespace qsc;

void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, int& var) {
  // toml11's find_or does not seem to give informative error messages
  // like it's "find" method does, hence my little function here.
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

/** Handle floats
 */
void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, qscfloat& var) {
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
void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, std::string& var) {
  if (indata.contains(varname)) {
    var = toml::find<std::string>(indata, varname);
  }
  varlist.push_back(varname);
}

/** Handle Vectors
 */
void toml_read(std::vector<std::string>& varlist, toml::value indata, std::string varname, Vector& var) {
  if (indata.contains(varname)) {
    auto indata_vector = toml::find<std::vector<double>>(indata, varname);
    // Convert the std::vector<double> to a Vector:
    var.resize(indata_vector.size(), 0.0);
    for (int j = 0; j < indata_vector.size(); j++) var[j] = (qscfloat)indata_vector[j];
  }
  varlist.push_back(varname);
}

/* Expand a Vector to a longer size, padding with zeros.
 */
void pad_vector(Vector& v, std::size_t newsize) {
  assert(v.size() <= newsize);
  Vector vcopy(v);

  v.resize(newsize, 0.0);
  for (int j = 0; j < vcopy.size(); j++) v[j] = vcopy[j];
}

/** Read in a configuration input file
 */
void Qsc::input(std::string filename) {
  if (verbose > 0) std::cout << "About to try loading input file " << filename << std::endl;
  
  auto indata = toml::parse(filename);

  std::vector<std::string> varlist;
  //auto varlist = std::vector<std::string>({"R0c", "R0s", "Z0c", "Z0s"});
  
  toml_read(varlist, indata, "nphi", nphi);
  toml_read(varlist, indata, "nfp", nfp);
  toml_read(varlist, indata, "eta_bar", eta_bar);
  toml_read(varlist, indata, "spsi", spsi);
  toml_read(varlist, indata, "sG", sG);
  toml_read(varlist, indata, "B0", B0);
  toml_read(varlist, indata, "I2", I2);
  toml_read(varlist, indata, "sigma0", sigma0);
  toml_read(varlist, indata, "B2s", B2s);
  toml_read(varlist, indata, "B2c", B2c);
  toml_read(varlist, indata, "p2", p2);
  toml_read(varlist, indata, "max_newton_iterations", max_newton_iterations);
  toml_read(varlist, indata, "max_linesearch_iterations", max_linesearch_iterations);
  toml_read(varlist, indata, "newton_tolerance", newton_tolerance);
  toml_read(varlist, indata, "verbose", verbose);
  toml_read(varlist, indata, "order_r_option", order_r_option);
  toml_read(varlist, indata, "R0c", R0c);
  toml_read(varlist, indata, "R0s", R0s);
  toml_read(varlist, indata, "Z0c", Z0c);
  toml_read(varlist, indata, "Z0s", Z0s);

  std::cout << "varlist:";
  int j;
  for (j = 0; j < varlist.size(); j++) std::cout << " " << varlist[j];
  std::cout << std::endl;

  // Check to see if there are any unused keys in the input file:
  int found_match;
  for (auto item : indata.as_table()) {
    // I believe "item" has type std::pair, representing a key-value pair.
    auto key = item.first;
    std::cout << "Found key: " << key << std::endl;
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

  // Now, make {R0c, R0s, Z0c, Z0s} have the same size.
  //auto sizes = std::valarray<int>({R0c.size(), R0s.size(), Z0c.size(), Z0s.size()});
  //int newsize = sizes.max();
  std::size_t newsize = 0;
  newsize = std::max(newsize, R0c.size());
  newsize = std::max(newsize, R0s.size());
  newsize = std::max(newsize, Z0c.size());
  newsize = std::max(newsize, Z0s.size());
  pad_vector(R0c, newsize);
  pad_vector(R0s, newsize);
  pad_vector(Z0c, newsize);
  pad_vector(Z0s, newsize);

  std::cout << "Arrays after reading input file:" << std::endl;
  std::cout << "R0c: " << R0c << std::endl;
  std::cout << "R0s: " << R0s << std::endl;
  std::cout << "Z0c: " << Z0c << std::endl;
  std::cout << "Z0s: " << Z0s << std::endl;
}
