#include <stdexcept>
#include <iostream>
#include <chrono>
#include <mpi.h>
#include "toml.hpp"
#include "qsc.hpp"
#include "scan.hpp"
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
  std::string infile(argv[1]);

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

  auto start = std::chrono::steady_clock::now();

  // Read general_option
  std::cout << "About to try loading input file " << infile << std::endl;
  auto indata = toml::parse(infile);
  std::vector<std::string> varlist;
  std::string general_option = GENERAL_OPTION_SINGLE;
  toml_read(varlist, indata, "general_option", general_option);

  if (general_option.compare(GENERAL_OPTION_SINGLE) == 0) {
    qsc::Qsc q;
    q.run(infile);

  } else if (general_option.compare(GENERAL_OPTION_RANDOM) == 0) {
    MPI_Init(NULL, NULL);
    
    qsc::Scan scan;
    scan.run(infile);

    MPI_Finalize();
    
  } else {
    throw std::runtime_error("Unrecognized setting for general_option");
  }
  
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Total time: "
            << elapsed.count() << " seconds" << std::endl;
  
  std::cout << "Good bye." << std::endl;
  
  return 0;
}
