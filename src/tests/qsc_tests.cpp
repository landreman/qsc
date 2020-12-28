#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#include <iostream>
#include <mpi.h>
#include "qsc.hpp"

// See https://github.com/onqtam/doctest/blob/master/doc/markdown/main.md
// main() taken from https://stackoverflow.com/questions/58289895/is-it-possible-to-use-catch2-for-testing-an-mpi-code 
int main( int argc, char* argv[] ) {
  if (qsc::single) {
    std::cout << "Running unit tests with SINGLE precision." << std::endl;
    // Make sure all output is written with full precision
    std::cout.precision(8);
    std::cerr.precision(8);
  } else {
    std::cout << "Running unit tests with DOUBLE precision." << std::endl;
    // Make sure all output is written with full precision
    std::cout.precision(15);
    std::cerr.precision(15);
  }
  
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  MPI_Init(&argc, &argv);
  int result = context.run();
  // int result = Catch::Session().run( argc, argv );
  MPI_Finalize();
  return result;
}
