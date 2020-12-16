#include <iostream>
#include "qsc.hpp"

int main(int argc, char* argv[]) {
  qsc::hw();

  qsc::Matrix m = qsc::differentiation_matrix(4, 0, 2 * qsc::pi);
  std::cout << m << std::endl;

  m = 4.3;
  std::cout << m(1, 1) << std::endl;
  std::cout << m << std::endl;
  
  return 0;
}
