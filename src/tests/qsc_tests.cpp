#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#include <mpi.h>

// See https://github.com/onqtam/doctest/blob/master/doc/markdown/main.md
// main() taken from https://stackoverflow.com/questions/58289895/is-it-possible-to-use-catch2-for-testing-an-mpi-code 
int main( int argc, char* argv[] ) {
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  MPI_Init(&argc, &argv);
  int result = context.run();
  // int result = Catch::Session().run( argc, argv );
  MPI_Finalize();
  return result;
}
