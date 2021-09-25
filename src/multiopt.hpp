#ifndef QSC_MULTIOPT_H
#define QSC_MULTIOPT_H

#include <vector>
#include "opt.hpp"

namespace qsc {

  class MultiOpt {
  private:    
    void defaults();
    
  public:
    std::vector<Opt> opts;
    int verbose;
    std::valarray<int> nphi;
    
    MultiOpt();
    void run(std::string);
    void input(std::string);
    void optimize();
    void write_netcdf();
  };
}

#endif

