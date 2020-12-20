#include <iostream>
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

  if (qsc::single) {
    std::cout << "Using SINGLE precision." << std::endl;
  } else {
    std::cout << "Using DOUBLE precision." << std::endl;
  }
  
  qsc::Qsc q;
    
  q.input(argv[1]);

  q.calculate();
  
  std::cout << "Good bye." << std::endl;
  
  return 0;
}
