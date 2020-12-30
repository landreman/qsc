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

