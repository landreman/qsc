#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include "qsc.hpp"

void qsc::hw() {
  std::cout << "Hello world!!!" << std::endl;

  int ierr;
  ierr = MPI_Init(NULL, NULL);
  if (ierr != 0) {
    throw std::runtime_error("Error in MPI_Init.");
  } else {
    std::cout << "MPI successfully initialized." << std::endl;
  }

  MPI_Finalize();
}

int qsc::return5() {
  return 5;
}
