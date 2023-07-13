# qsc - Quasisymmetric Stellarator Construction

This code generates quasisymmetric stellarator configurations using a
near-axis expansion, as detailed in [Landreman \& Sengupta,
"Constructing stellarators with quasisymmetry to high order",
J. Plasma Phys. 85, 815850601
(2019)](http://doi.org/10.1017/S0022377819000783).  This repository
contains a C++ implementation of the method. For a python
implementation, see [pyQSC](https://github.com/landreman/pyQSC). For a
fortran implementation, see
[here](https://github.com/landreman/quasisymmetry). The C++ version
here was used for the paper [Landreman, "Mapping the space of
quasisymmetric stellarators using optimized near-axis expansion",
J. Plasma Phys. 88, 905880616
(2022)](http://doi.org/10.1017/S0022377822001258). This C++ version is
best if you are doing large parameter scans or optimization, where
speed is a priority. For use cases in which you can afford to spend
100 ms rather than 1ms to compute an equilibrium, you should consider
[pyQSC](https://github.com/landreman/pyQSC).

This code provides a standalone executable `xqsc`, which reads input
files and writes output files. However a library `libqsc.a` is also
generated, and you can write your own driver programs that link to
this library.


## Requirements
* A C++ compiler (C++11 standard)
* CMake
* BLAS
* LAPACK
* MPI
* NetCDF (The C++ interface to NetCDF is not required, just the standard C interface.)


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
and the driver program `xqsc` will be present in the `bin/` directory.

If you want to clear all CMake files to do a clean configure and build, you
can run the `./clean_cmake` script.

### Single precision

By default the code uses double precision for all calculations, but it
is also possible to build the code using single precision. To do this,
run
~~~~
cmake -DSINGLE=ON .
~~~~
(Alternative flags like `-DSINGLE=1` also work.) Then run `make` as
before.  The library and driver that are then compiled will be named
`libqsc_single.a` and `xqsc_single`. I have not seen an advantage in
practice to using single precision, so you almost certainly want to
use the standard double-precision build instead.


## Testing

After CMake has been run, the unit tests can be built and run by typing
~~~~
make test
~~~~
It is also possible to type `make unitTests` to build the tests without running them.
The unit test executable is named `unitTests` and is placed in the project's `tests/` directory.


## Running the code

To run the code, type
~~~~
<path-to-executable>/xqsc qsc_in.<extension>
~~~~
where `<path-to-executable>` is usually `bin`. An output file will be
generated, with the name `qsc_out.<extension>.nc`. The output files
are in NetCDF format. You can browse the contents of the output files
using `ncdump qsc_out.<extension>.nc | less`.


### Input files

Input files should have filenames of the form
`qsc_in.<extension>`. The input files use TOML format. Many example
input files can be found in the `examples` directory. Please consult
these examples for the typical variables that should be specified for
each type of calculation. Every input file has a top-level option
`general_option` which controls the type of calculation that is run.

### Settings for general_option

#### `general_option = "single"`

For this setting, `xqsc` will solve the near-axis equations for a
single configuration. The input parameters for this configuration are
specified in the `qsc` section of the TOML input file (i.e. under
`[qsc]`.)

A calculation with `general_option = "single"` typically takes about 1
ms.

MPI is not used for this setting.

#### `general_option = "opt"`

For this setting, a single optimization will be run. The initial
configuration is specified in the `[qsc]` section of the input
file. The additional settings that determine the objective function
and parameter space are specified in the `[opt]` section of the input
file. One set of these parameters determines which variables are fixed
or varied in the optimization, e.g. `vary_eta_bar` and
`vary_R0c`. Another set of parameters determines which terms are
included in the objective function, and if so, how strongly they are
weighted. Examples of these variables are `weight_grad_B` and
`weight_B20`. For any weights that are ≤ 0, the corresponding term is
not included in the objective.

Optimizations typically use Fourier refinement, in which the number of
Fourier modes in the axis shape is increased in stages. The number of
steps of refinement is specified by the `fourier_refine` variable, an
integer. If `fourier_refine = 2`, then optimization will be run with
N, N + 1, and N + 2 sine and cosine modes in R and Z, where N is the
number of sine/cosine modes in R/Z in the `[qsc]` list. The phi
resolution for each stage can be specified in the `nphi` list.

A calculation with `general_option = "opt"` typically takes on the
order of a second or a few seconds, depending on resolution and the
number of Fourier refinement steps.

MPI is not used for this setting.

#### `general_option = "multiopt"`

For this setting, several stages of optimization are applied to one
configuration. Each stage generally has a different objective
function. In typical usage, you might have two multiopt stages. In the
first stage, the objective includes just the ∇B scale length and B20,
whereas the second stage also includes other objectives such as the
∇∇B scale length.

For multiopt runs, a `[multiopt]` section of the input file is
required. In this section, a variable `nopts` specifies how many opt
stages are to be performed. Each opt stage is exactly like a
calculation with `general_option = "opt"`, except the `[opt]` section
in the input file is replaced with `[opt0]`, `[opt1]`, etc for the
different stages.

As with `general_option = "opt"`, the initial condition for the
optimization is determined by the `[qsc]` list in the input file.

A calculation with `general_option = "multiopt"` typically takes on the
order of a second or a few seconds, depending on resolution and the
number of Fourier refinement steps.

MPI is not used for this setting.

#### `general_option = "scan"`

For this setting, a large number of configurations are computed. No
optimization is performed. Input parameters are chosen randomly. MPI
is used, and each process performs calculations independently of the
others. Solutions are discarded if they do not meet certain criteria
("filters"), e.g. if the minimum aspect ratio is too large. In this
type of scan, often the yield is low, in the sense that most
configurations get filtered out.

MPI is used for this setting, so the executable should be launched
appropriately, e.g. `mpiexec -n 4 xqsc qsc_in.random_scan_small` or
`srun xqsc qsc_in.random_scan_small`. Running with just one MPI
process will not cause an error, but you probably want to take
advantage of MPI for large scans.

For this setting, the scan settings are defined in a `[scan]` section
of the input file. Here, the minimum and maximum values for each input
parameter are defined. You can also specify if each input parameter is
distributed uniformly or logarithmically. The filters are also defined
in this section. Examples of variables defining the filters are
`min_iota_to_keep`, `max_B20_variation_to_keep`, and
`min_r_singularity_to_keep`. The scan will run for a time specified by
the `max_seconds` variable.

For this setting, the `[qsc]` list is still required. This list
specifies all the input parameters for each configuration other than
the ones that are set randomly during the scan.

#### `general_option = "multiopt_scan"`

For this setting, multiple "multiopt" optimizations are run. Each
optimization typically has different weights or target values for the
various terms in the objective function. Solutions are discarded if
they do not meet certain criteria ("filters"), e.g. if the minimum
aspect ratio is too large. Compared to `general_option = "scan"`, the
scans in `general_option = "multiopt_scan"` typically have higher
yield, i.e. more configurations pass the filters. This is because each
configuration is the result of an optimization rather than the result
of random input parameters. However each configuration requires an
optimization, so more time is needed for each candidate configuration.

This `general_option = "multiopt_scan"` setting is the one used for
the large parameter scans described in [Landreman, "Mapping the space
of quasisymmetric stellarators using optimized near-axis expansion",
J. Plasma Phys. 88, 905880616
(2022)](http://doi.org/10.1017/S0022377822001258).

MPI is used for this setting, so the executable should be launched
appropriately, e.g. `mpiexec -n 4 xqsc
qsc_in.multiopt_scan_QH_nfp4_small_keep_all_1` or `srun xqsc
qsc_in.multiopt_scan_QH_nfp4_small_keep_all_1`. Running with just one
MPI process will cause an error - you must use at least two processes.
This is because process 0 is reserved for coordinating the scan, and
at least one other process is needed to actually do computation.

A `[multiopt_scan]` section of the input file contains the parameters
that define the scan. The key variables are `params`, `params_min`,
`params_max`, `params_n`, `params_log`, and `params_stage`. These
variables are each lists of the same length. Here, `params` is a list
of strings, indicating which parameters to vary. A tensor-product scan
is performed, in which the values considered for each variable are
independent of the values of the other variables. The maximum and
minimum values for each variable are determined by `params_min` and
`params_max`. The number of values for each variable is set by
`params_n`. The list `params_log` contains booleans, determining if
the values for each variable are spaced logarithmically (true) or
linearly (false). The variable `params_stage` determines which stage
of the multiopt the corresponding variable applies to; a value of -1
indicates all stages.

The `[multiopt_scan]` section of the input file also contains the
filter variables, e.g. `max_B20_variation_to_keep`.  The input file
should contain all the lists for `general_option = "multiopt"`, such
as `[multiopt]`, `[opt0]`, and `[opt1]`. Also the `[qsc]` list is
still required, specifying the initial conditions for each
optimization.

Note that the wallclock time for a multiopt_scan depends on how many
values you choose to consider in the scan.


## Interaction with pyQSC

For a `qsc_out.<extension>.nc` file generated by this C++ code, you
can load the result into [pyQSC](https://github.com/landreman/pyQSC)
using the `from_cxx()` method, for plotting and other analysis:
~~~~py
from qsc import Qsc

config = Qsc.from_cxx("qsc_out.my_configuration.nc")
~~~~

This C++ repository provides the `qscPlot` script in the `bin`
directory, which makes use of pyQSC for plotting configurations
generated in C++.