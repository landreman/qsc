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

If CMake has any trouble finding the required depencies, or if you
want to tell it to use a different version of a compiler or library
than the version it detected, you have two options. One option is to
provide a flag when you call CMake, e.g. `cmake
-DCMAKE_CXX_COMPILER=mpiicpc .` The second option is to edit the file
`cmake/knownHosts.cmake` to add a case for your machine. In this file,
CMake attempts to detect if it is being run on a known machine based
on the environment variables, and if it recognizes the machine,
appropriate settings for that machine are set.

Once CMake has completed successfully, then run
~~~~
make
~~~~
If the build is successful, the `qsc` library will be built in the project's `lib/` directory,
and `qsc_driver` will be present in the `bin/` directory.

If you want to clear all CMake files to do a clean configure and build, you
can run the `./clean_cmake` script.

### Single precision

By default the code uses double precision for all calculations, but it
is also possible to build the code using single precision. To do this,
run
~~~~
cmake -DSINGLE=1 .
~~~~
(Alternative flags like `-DSINGLE=ON` also work.) Then run `make` as before.
The library and driver that are then compiled will be named `libqsc_single.a` and `qsc_driver_single`.


## Testing

After CMake has been run, the unit tests can be built and run by typing
~~~~
make test
~~~~
It is also possible to type `make unitTests` to build the tests without running them.
The unit test executable is named `unitTests` and is placed in the project's `tests/` directory.

