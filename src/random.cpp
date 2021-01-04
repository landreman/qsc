#include <cmath>
#include <chrono>
#include <random>
#include <stdexcept>
#include "random.hpp"

using namespace qsc;

/**
 * Policy on signs:
 * - Swapping min and max never has any effect.
 * - For "linear", either sign is allowed for min and max.
 * - For "log", min and max must have the same sign, either + or -. 
 *   If min and max are negative, values returned by get() will be negative.
 *   It is also valid to have min=max=0, in which case get() always returns 0.
 * - For "2 sided log", min and max must both be positive, to avoid confusion.
 */
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
  logmin = log(std::abs(min));
  logmax = log(std::abs(max));
  last = 0.0;

  sign = 1;
  if (max < 0) sign = -1;

  unsigned my_seed = std::chrono::steady_clock::now().time_since_epoch().count();
  generator.seed(my_seed);
  
  // Convert string to int for speed later.
  if (distrib_in.compare(RANDOM_OPTION_LINEAR) == 0) {
    distrib = RANDOM_INT_OPTION_LINEAR;
    
  } else if (distrib_in.compare(RANDOM_OPTION_LOG) == 0) {
    distrib = RANDOM_INT_OPTION_LOG;
    if (max * min < 0) throw std::runtime_error("Error initializing Random with log distribution: max and min must have the same sign");
    if ((max != min) && (max == 0 || min == 0))
      throw std::runtime_error("Error initializing Random with log distribution: max and min must both be nonzero, unless they are both 0");
    
  } else if (distrib_in.compare(RANDOM_OPTION_2_SIDED_LOG) == 0) {
    distrib = RANDOM_INT_OPTION_2_SIDED_LOG;
    if (min == 0 && max == 0) {
      // Ok, get() will return 0
    } else if (min <= 0 || max <= 0) {
      throw std::runtime_error("Error initializing Random with 2-sided log distribution: max and min must both be positive, or both 0");
    }
    
  } else {
    throw std::runtime_error("Unrecognized random distribution type");
  }
}

qscfloat Random::get() {
  // Explicitly handle the case of min==max, to avoid evaluating log(0).
  if (max == min) return max;
  //if (std::abs(max - min) < 1.0e-30) return max;
  
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
    val = sign * exp(logmin + (logmax - logmin) * rand01);
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

/** For deterministic=true, set the position in the "random" number
 * sequence to the 0-based position n. For deterministic=false this
 * function has no effect.
 */
void Random::set_to_nth(int n) {
  // This next version introduces roundoff-scale errors that are
  //sufficient to affect single-precision tests.
  // last = std::fmod(n * gamma, (qscfloat)1.0);
  
  last = 0.0;
  for (int j = 0; j < n; j++) get();
}
