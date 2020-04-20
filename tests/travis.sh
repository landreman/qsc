#!/bin/bash

#set -ex
set -e
# In the above line, "set -x" causes this script to exit as soon as any line fails. "set -e" causes each line of this script to be printed (with a + in front) before it is executed, so if a step fails, you can see from the travis log what command failed.

echo Hello from travis.sh


pwd
ls
env
which mpicc
which mpiexec
which mpirun

cmake --version

cmake .

make

bin/qsc_driver

mpiexec -n 2 bin/qsc_driver