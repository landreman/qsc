#include <stdexcept>
#include "doctest.h"
#include "random.hpp"

using doctest::Approx;

TEST_CASE("Test Random class for linear distribution") {
  qsc::qscfloat min = -3.0, max = -2.0;
  int N = 30000;
  int n_intervals = 10;
  int histogram[n_intervals];
  int j, bin;
  qsc::qscfloat val;

  for (int deterministic = 0; deterministic < 2; deterministic++) {
    qsc::Random r((bool)deterministic, "linear", min, max);

    for (j = 0; j < n_intervals; j++) histogram[j] = 0;
    
    for (j = 0; j < N; j++) {
      val = r.get();
      bin = (int)floor((val - min) / (max - min) * n_intervals);
      histogram[bin]++;
      //std::cout << val << " " << bin << " | ";
    }
    //std::cout << std::endl;
    
    std::cout << "Final histogram: ";
    for (j = 0; j < n_intervals; j++) std::cout << histogram[j] << " ";
    std::cout << std::endl;
    
    // There should be about N / n_intervals hits in each bin:
    for (j = 0; j < n_intervals; j++) {
      CHECK(histogram[j] > (0.9 * N) / n_intervals);
      CHECK(histogram[j] < (1.1 * N) / n_intervals);
    }
  }
}

TEST_CASE("Test Random class for log distribution") {
  qsc::qscfloat min = 1.0e+1, max = 1.0e+9;
  int N = 30000;
  int n_intervals = 8;
  int histogram[n_intervals];
  int j, bin;
  qsc::qscfloat val;

  for (int deterministic = 0; deterministic < 2; deterministic++) {
    qsc::Random r((bool)deterministic, "log", min, max);

    for (j = 0; j < n_intervals; j++) histogram[j] = 0;
    
    for (j = 0; j < N; j++) {
      val = r.get();
      //std::cout << val << " ";
      if (val >= 1.0e+1 && val < 1.0e+2) histogram[0]++;
      if (val >= 1.0e+2 && val < 1.0e+3) histogram[1]++;
      if (val >= 1.0e+3 && val < 1.0e+4) histogram[2]++;
      if (val >= 1.0e+4 && val < 1.0e+5) histogram[3]++;
      if (val >= 1.0e+5 && val < 1.0e+6) histogram[4]++;
      if (val >= 1.0e+6 && val < 1.0e+7) histogram[5]++;
      if (val >= 1.0e+7 && val < 1.0e+8) histogram[6]++;
      if (val >= 1.0e+8 && val < 1.0e+9) histogram[7]++;
    }
    //std::cout << std::endl;
    
    std::cout << "Final histogram: ";
    for (j = 0; j < n_intervals; j++) std::cout << histogram[j] << " ";
    std::cout << std::endl;
    
    // There should be about N / n_intervals hits in each bin:
    for (j = 0; j < n_intervals; j++) {
      CHECK(histogram[j] > (0.9 * N) / n_intervals);
      CHECK(histogram[j] < (1.1 * N) / n_intervals);
    }
  }
}


TEST_CASE("Test Random class for 2-sided log distribution") {
  qsc::qscfloat min = 1.0e+1, max = 1.0e+5;
  int N = 30000;
  int n_intervals = 8;
  int histogram[n_intervals];
  int j, bin;
  qsc::qscfloat val;

  for (int deterministic = 0; deterministic < 2; deterministic++) {
    qsc::Random r((bool)deterministic, "2 sided log", min, max);

    for (j = 0; j < n_intervals; j++) histogram[j] = 0;
    
    for (j = 0; j < N; j++) {
      val = r.get();
      // std::cout << val << " ";
      if (val <= -1.0e+1 && val > -1.0e+2) histogram[0]++;
      if (val <= -1.0e+2 && val > -1.0e+3) histogram[1]++;
      if (val <= -1.0e+3 && val > -1.0e+4) histogram[2]++;
      if (val <= -1.0e+4 && val > -1.0e+5) histogram[3]++;
      if (val >= 1.0e+1 && val < 1.0e+2) histogram[4]++;
      if (val >= 1.0e+2 && val < 1.0e+3) histogram[5]++;
      if (val >= 1.0e+3 && val < 1.0e+4) histogram[6]++;
      if (val >= 1.0e+4 && val < 1.0e+5) histogram[7]++;
    }
    // std::cout << std::endl;
    
    std::cout << "Final histogram: ";
    for (j = 0; j < n_intervals; j++) std::cout << histogram[j] << " ";
    std::cout << std::endl;
    
    // There should be about N / n_intervals hits in each bin:
    for (j = 0; j < n_intervals; j++) {
      CHECK(histogram[j] > (0.9 * N) / n_intervals);
      CHECK(histogram[j] < (1.1 * N) / n_intervals);
    }
  }
}

TEST_CASE("Test Random::set_to_nth") {
  qsc::qscfloat min = -3.0, max = -2.0;
  bool deterministic = true;
  qsc::Random r1(deterministic, "linear", min, max);
  qsc::Random r2(deterministic, "linear", min, max);
  int pos = 17;
  for (int j = 0; j < pos; j++) r1.get();
  r2.set_to_nth(pos);
  CHECK(Approx(r1.get()) == r2.get());
}

TEST_CASE("Random: verify results make sense if min == max") {
  /* The tricky case is a 2-sided log. One reasonable policy for this
     distribution if min==max would be to return either min or -min
     with equal probability. Instead, I have chosen the policy that
     you return only min, not -min.
   */
  qsc::qscfloat min, max;
  bool deterministic;
  int k, j;

  // Try negative, zero, and positive values
  for (j = 0; j < 3; j++) {
    switch (j) {
    case 0:
      min = -3.14;
      break;
    case 1:
      min = 0;
      break;
    case 2:
      min = 3.14;
      break;
    default:
      throw std::runtime_error("Should not get here");
    }
    max = min;
    for (int j_deterministic = 0; j_deterministic < 2; j_deterministic++) {
      deterministic = (bool) j_deterministic;
      CAPTURE(deterministic);
      
      qsc::Random r1(deterministic, "linear", min, max);
      for (k = 0; k < 5; k++) CHECK(Approx(r1.get()) == min);
      
      qsc::Random r2(deterministic, "log", min, max);
      for (k = 0; k < 5; k++) CHECK(Approx(r2.get()) == min);

      if (min < 0) {
	CHECK_THROWS((new qsc::Random(deterministic, "2 sided log", min, max)));
      } else {
	// Both 0 and positive min/max should work:
	qsc::Random r3(deterministic, "2 sided log", min, max);
	for (k = 0; k < 5; k++) CHECK(Approx(r3.get()) == min);
      }
    }
  }
}

TEST_CASE("Random: verify results make sense for negative min and/or max") {
  qsc::qscfloat min, max, val;
  int j;
  bool deterministic, complicated_bool;

  for (int j_deterministic = 0; j_deterministic < 2; j_deterministic++) {
    deterministic = (bool) j_deterministic;
    CAPTURE(deterministic);
    
    // Log, both negative
    min = -2;
    max = -1;
    qsc::Random r1(deterministic, "log", min, max);
    for (j = 0; j < 10; j++) {
      val = r1.get();
      CHECK(val >= -2);
      CHECK(val <= -1);
    }
    // Same but min and max swapped:
    min = -1;
    max = -2;
    qsc::Random r2(deterministic, "log", min, max);
    for (j = 0; j < 10; j++) {
      val = r2.get();
      CHECK(val >= -2);
      CHECK(val <= -1);
    }

    // Log with one positive and one negative should throw an exception.
    CHECK_THROWS((new qsc::Random(deterministic, "log", -1.0, 2.0)));
    // It is an error if either min or max is 0, but not both
    CHECK_THROWS((new qsc::Random(deterministic, "log", 0, 2.0)));
    CHECK_THROWS((new qsc::Random(deterministic, "log", 0, -2.0)));
    CHECK_THROWS((new qsc::Random(deterministic, "log", 2.0, 0)));
    CHECK_THROWS((new qsc::Random(deterministic, "log", -2.0, 0)));
    
    // Linear, both negative
    min = -2;
    max = -1;
    qsc::Random r3(deterministic, "linear", min, max);
    for (j = 0; j < 10; j++) {
      val = r3.get();
      CHECK(val >= -2);
      CHECK(val <= -1);
    }
    // Same but min and max swapped:
    min = -1;
    max = -2;
    qsc::Random r4(deterministic, "linear", min, max);
    for (j = 0; j < 10; j++) {
      val = r4.get();
      CHECK(val >= -2);
      CHECK(val <= -1);
    }

    // Linear, one of each sign
    min = -2;
    max = 0.1;
    qsc::Random r5(deterministic, "linear", min, max);
    for (j = 0; j < 10; j++) {
      val = r5.get();
      CHECK(val >= -2);
      CHECK(val <= 0.1);
    }
    // Same but min and max swapped:
    min = 0.1;
    max = -2;
    qsc::Random r6(deterministic, "linear", min, max);
    for (j = 0; j < 10; j++) {
      val = r6.get();
      CHECK(val >= -2);
      CHECK(val <= 0.1);
    }
    
    // 2-sided log:
    // One positive and one negative bound:
    CHECK_THROWS((new qsc::Random(deterministic, "2 sided log", -1.0, 2.0)));
    CHECK_THROWS((new qsc::Random(deterministic, "2 sided log", 1.0, -2.0)));
    // 2 negative unequal bounds:
    CHECK_THROWS((new qsc::Random(deterministic, "2 sided log", -1.0, -2.0)));
    CHECK_THROWS((new qsc::Random(deterministic, "2 sided log", -2.0, -1.0)));
    // Equal negative bounds:
    CHECK_THROWS((new qsc::Random(deterministic, "2 sided log", -2.0, -2.0)));
    
  }
}
