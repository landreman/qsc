#include <exception>
#include "qsc.hpp"
#include "toml.hpp"

using namespace qsc;

void toml_read(std::string filename, toml::value indata, std::string varname, int& var) {
  // toml11's find_or does not seem to give informative error messages
  // like it's "find" method does, hence my little function here.
  if (indata.contains(varname)) {
    var = toml::find<int>(indata, varname);
  }
  /*
  try {
    var = toml::find_or(indata, varname, var);
  } catch (std::exception) {
    std::cerr << "Unable to read variable " << varname << " from input file " << filename << std::endl;
    throw;
  }
  */
}

void toml_read(std::string filename, toml::value indata, std::string varname, qscfloat& var) {
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
}

void Qsc::input(std::string filename) {
  if (verbose > 0) std::cout << "About to try loading input file " << filename << std::endl;
  
  auto indata = toml::parse(filename);

  toml_read(filename, indata, "nphi", nphi);
  toml_read(filename, indata, "nfp", nfp);
  toml_read(filename, indata, "eta_bar", eta_bar);
  toml_read(filename, indata, "spsi", spsi);
  toml_read(filename, indata, "sG", sG);
  toml_read(filename, indata, "B0", B0);
  toml_read(filename, indata, "I2", I2);
  toml_read(filename, indata, "sigma0", sigma0);
  toml_read(filename, indata, "B2s", B2s);
  toml_read(filename, indata, "B2c", B2c);
  toml_read(filename, indata, "p2", p2);
  toml_read(filename, indata, "max_newton_iterations", max_newton_iterations);
  toml_read(filename, indata, "max_linesearch_iterations", max_linesearch_iterations);
  toml_read(filename, indata, "newton_tolerance", newton_tolerance);
  toml_read(filename, indata, "verbose", verbose);
}
