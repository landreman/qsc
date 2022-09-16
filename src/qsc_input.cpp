#include <chrono>
#include "qsc.hpp"
#include "toml.hpp"
#include "toml_util.hpp"

using namespace qsc;

/** Read in a configuration input file
 */
void qsc::Qsc::input(std::string filename) {
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) start = std::chrono::steady_clock::now();
  
  auto toml_file = toml::parse(filename);
  auto indata = toml::find(toml_file, "qsc");
  
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
  toml_read(varlist, indata, "fc", fc);
  toml_read(varlist, indata, "fs", fs);
  toml_read(varlist, indata, "angle_shift", angle_shift);

  toml_unused(varlist, indata);
  
  // Now, make {R0c, R0s, Z0c, Z0s} have the same size.
  std::size_t newsize = 0;
  newsize = std::max(newsize, R0c.size());
  newsize = std::max(newsize, R0s.size());
  newsize = std::max(newsize, Z0c.size());
  newsize = std::max(newsize, Z0s.size());
  newsize = std::max(newsize, fc.size());
  newsize = std::max(newsize, fs.size());
  pad_vector(R0c, newsize);
  pad_vector(R0s, newsize);
  pad_vector(Z0c, newsize);
  pad_vector(Z0s, newsize);
  pad_vector(fc, newsize);
  pad_vector(fs, newsize);

  std::cout << "Fourier amplitudes for the magnetic axis shape:" << std::endl;
  std::cout << "R0c: " << R0c << std::endl;
  std::cout << "R0s: " << R0s << std::endl;
  std::cout << "Z0c: " << Z0c << std::endl;
  std::cout << "Z0s: " << Z0s << std::endl;
  std::cout << "fc:  " <<  fc << std::endl;
  std::cout << "fs:  " <<  fs << std::endl;

  if (verbose > 0) {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for reading input: "
	      << elapsed.count() << " seconds" << std::endl;
  }
  
}
