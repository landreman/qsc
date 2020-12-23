#include <string>
#include <stdexcept>
#include "doctest.h"
#include "qsc.hpp"

using namespace qsc;
using doctest::Approx;

TEST_CASE("outfile") {
  CHECK(qsc::outfile("qsc_in.ext") == "qsc_out.ext.nc");
  CHECK(qsc::outfile("/path/to/infile/qsc_in.ext") == "/path/to/infile/qsc_out.ext.nc");
  CHECK(qsc::outfile("/qsc_in.ext") == "/qsc_out.ext.nc");

  /*
    // doctest CHECK_THROWS does not seem to be working.
  std::string mystr = "foo";
  CHECK_THROWS_AS(qsc::outfile(mystr), std::runtime_error);
  CHECK_THROWS(qsc::outfile("foobarbar"));
  CHECK_THROWS(qsc::outfile("foobar"));
  CHECK_THROWS(qsc::outfile("qsc_in")); // Missing period at end
  CHECK_THROWS(qsc::outfile("/path/qsc_in.foo/bar"));
  CHECK_THROWS(qsc::outfile("/path/qsc_in.foo/"));
  */
}

