#include <iostream>
#include "qsc.hpp"

int main(int argc, char* argv[]) {
  qsc::hw();

  qsc::Vec m = qsc::differentiation_matrix(3, 0, 2 * qsc::pi);
  
  return 0;
}
