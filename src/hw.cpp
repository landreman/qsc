#include <iostream>
#include <mpi.h>
#include "qsc.hpp"

void qsc::hw() {
  std::cout << "Hello world!!!" << std::endl;

  int ierr;
  ierr = MPI_Init(NULL, NULL);
  if (ierr != 0) {
    std::cerr << "Error in MPI_Init." << std::endl;
    exit(1);
  } else {
    std::cout << "MPI successfully initialized." << std::endl;
  }

  MPI_Finalize();
}

int qsc::return5() {
  return 5;
}
