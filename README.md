# qsc
Quasisymmetric Stellarator Construction

## Requirements
* A C++ compiler
* CMake
* BLAS
* LAPACK
* MPI

## Building the code

From the root project directory, run
~~~~
cmake .
~~~~
Then run
~~~~
make
~~~~
If the build is successful, the `qsc` library will be built in the project's `lib/` directory,
and `qsc_driver` will be present in the `bin/` directory.

## Testing

After CMake has been run, the unit tests can be built and run by typing
~~~~
make test
~~~~
It is also possible to type `make unitTests` to build the tests without running them.
The unit test executable is named `unitTests` and is placed in the project's `tests/` directory.