#!/bin/bash

echo
echo "All of Catch2's code is directly included in this respository, so there is usually no need to run this script,"
echo "unless you want to update Catch2 to the latest version. In this case, first delete the existing catch2 directory, then run this script."
echo

# For discussion of including the single-file version of Catch2 directly in repositories, see
# https://github.com/catchorg/Catch2/blob/master/docs/tutorial.md#where-to-put-it
# and
# https://levelofindirection.com/blog/unit-testing-in-cpp-and-objective-c-just-got-ridiculously-easier-still.html

set -ex
# In the above line, "set -e" causes this script to exit as soon as any line fails. "set -x" causes each line of this script to be printed (with a + in front) before it is executed, so if a step fails, you can see from the travis log what command failed.

echo Hello from install_catch2.sh

pwd

mkdir catch2
cd catch2
# The following command gets the latest version:
wget https://raw.githubusercontent.com/catchorg/Catch2/master/single_include/catch2/catch.hpp
# Or, you can use the following command to get a specific version:
#wget https://github.com/catchorg/Catch2/releases/download/v2.11.3/catch.hpp
cd ..

echo Done installing Catch2
