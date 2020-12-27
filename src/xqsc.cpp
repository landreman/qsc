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

  if (q.at_least_order_r2) {
    auto grad_grad_B_alt = q.calculate_grad_grad_B_tensor_alt();
    for (int j2 = 0; j2 < 3; j2++) {
      for (int j3 = 0; j3 < 3; j3++) {
	for (int j4 = 0; j4 < 3; j4++) {
	  std::cout << "grad_grad_B(0," << j2 << "," << j3 << "," << j4 << ")="
		    << q.grad_grad_B_tensor(0, j2, j3, j4) << " or "
		    << grad_grad_B_alt(0, j2, j3, j4) << std::endl;
	}
      }
    }
  }

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
