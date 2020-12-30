#include <cmath>
#include <random>
#include <stdexcept>
#include "random.hpp"

using namespace qsc;

Random::Random(bool deterministic_in,
	       std::string distrib_in,
	       qscfloat min_in,
	       qscfloat max_in) :
  // Call constructor of member objects:
  uniform_distrib(0.0, 1.0)
{
  deterministic = deterministic_in;
  min = min_in;
  max = max_in;
  logmin = log(min);
  logmax = log(max);
  last = 0.0;

  // Convert string to int for speed later.
  if (distrib_in.compare(RANDOM_OPTION_LINEAR) == 0) {
    distrib = RANDOM_INT_OPTION_LINEAR;
  } else if (distrib_in.compare(RANDOM_OPTION_LOG) == 0) {
    distrib = RANDOM_INT_OPTION_LOG;
  } else if (distrib_in.compare(RANDOM_OPTION_2_SIDED_LOG) == 0) {
    distrib = RANDOM_INT_OPTION_2_SIDED_LOG;
  } else {
    throw std::runtime_error("Unrecognized random distribution type");
  }
}

qscfloat Random::get() {
  // First, get a random float in [0, 1]:
  qscfloat rand01, temp;
  if (deterministic) {
    last = std::fmod(last + gamma, (qscfloat)1.0);
    rand01 = last;
  } else {
    rand01 = uniform_distrib(generator);
  }

  // Map to the distribution we want
  qscfloat val = 0.0;
  switch (distrib) {
  case RANDOM_INT_OPTION_LINEAR:
    val = min + (max - min) * rand01;
    break;

  case RANDOM_INT_OPTION_LOG:
    val = exp(logmin + (logmax - logmin) * rand01);
    break;

  case RANDOM_INT_OPTION_2_SIDED_LOG:
    temp = std::abs((rand01 - 0.5) * 2);
    val = exp(logmin + (logmax - logmin) * temp);
    if (rand01 > 0.5) val = -val;
    break;

  default:
    throw std::runtime_error("Unrecognized random distribution type");
  }
  
  return val;
}
