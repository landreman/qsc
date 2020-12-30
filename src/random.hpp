#ifndef QSC_RANDOM_H
#define QSC_RANDOM_H

#include <random>
#include "vector_matrix.hpp"

namespace qsc {

  const std::string RANDOM_OPTION_LINEAR = "linear";
  const std::string RANDOM_OPTION_LOG = "log";
  const std::string RANDOM_OPTION_2_SIDED_LOG = "2 sided log";

  enum Distrib_type {
    RANDOM_INT_OPTION_LINEAR,
    RANDOM_INT_OPTION_LOG,
    RANDOM_INT_OPTION_2_SIDED_LOG};
  
  class Random {
  private:
    qscfloat min, max, logmin, logmax;
    bool deterministic;
    Distrib_type distrib;
    std::default_random_engine generator;
    std::uniform_real_distribution<qscfloat> uniform_distrib;
    const qscfloat gamma = (1.0 + sqrt(5.0)) / 2.0;
    qscfloat last;
    
  public:    
    Random(bool, std::string, qscfloat, qscfloat);
    qscfloat get();
  };
}

#endif

