#include <iostream>
#include <ctime>
#include <chrono>
#include "qsc.hpp"

int main(int argc, char* argv[]) {
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
  } else {
    std::cout << "Using DOUBLE precision." << std::endl;
  }

  std::time_t start_time, end_time;
  start_time = std::clock();
  auto start = std::chrono::steady_clock::now();
  
  qsc::Qsc q;
    
  q.input(directory_and_infile);

  q.calculate();

  q.write_netcdf(directory_and_outfile);
  
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
