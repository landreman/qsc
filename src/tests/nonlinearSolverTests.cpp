#include "doctest.h"
#include "qsc.hpp"

TEST_CASE("foo") {
  CHECK(qsc::return5() + 1 == 6);
}
