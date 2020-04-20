#include "catch.hpp"
#include "qsc.hpp"

TEST_CASE("foo","[bar]") {
  CHECK(qsc::return5() + 1 == 6);
}
