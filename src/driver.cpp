#include <stdexcept>
#include <iostream>
#include <ctime>
#include <chrono>
#include "toml.hpp"
#include "qsc.hpp"
#include "toml_util.hpp"

const std::string GENERAL_OPTION_SINGLE = "single";
const std::string GENERAL_OPTION_RANDOM = "random";

int qsc::driver(int argc, char* argv[]) {
  std::cout << "QSC: Quasisymmetric Stellarator Construction" << std::endl;

  if (argc != 2) {
    std::string exe;
    if (qsc::single) {
      exe = "xqsc_single";
    } else {
      exe = "xqsc";
    }
    std::cout << "Usage: " << exe << " qsc_in.<extension>" << std::endl;
    return 1;
  }
  std::string directory_and_infile(argv[1]);
  std::string directory_and_outfile = qsc::outfile(directory_and_infile);

  if (qsc::single) {
    std::cout << "Using SINGLE precision." << std::endl;
    // Make sure all output is written with full precision
    std::cout.precision(8);
    std::cerr.precision(8);
  } else {
    std::cout << "Using DOUBLE precision." << std::endl;
    // Make sure all output is written with full precision
    std::cout.precision(15);
    std::cerr.precision(15);
  }

  std::time_t start_time, end_time;
  start_time = std::clock();
  auto start = std::chrono::steady_clock::now();

  // Read general_option
  std::cout << "About to try loading input file " << directory_and_infile << std::endl;
  auto indata = toml::parse(directory_and_infile);
  std::vector<std::string> varlist;
  std::string general_option = GENERAL_OPTION_SINGLE;
  toml_read(varlist, indata, "general_option", general_option);

  if (general_option.compare(GENERAL_OPTION_SINGLE) == 0) {
    qsc::Qsc q;    
    q.input(indata);
    q.calculate();
    q.write_netcdf(directory_and_outfile);

  } else if (general_option.compare(GENERAL_OPTION_RANDOM) == 0) {
    std::cout << "Random scan goes here." << std::endl;
    
  } else {
    throw std::runtime_error("Unrecognized setting for general_option");
  }
  
  end_time = std::clock();
  auto end = std::chrono::steady_clock::now();

  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Total time from chrono:           "
            << elapsed.count() << " seconds" << std::endl;
  std::cout << "Total time from ctime (CPU time): "
            << double(end_time - start_time) / CLOCKS_PER_SEC
	    << " seconds" << std::endl;
  
  std::cout << "Good bye." << std::endl;
  
  return 0;
}
