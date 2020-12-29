#ifndef QSC_SCAN_H
#define QSC_SCAN_H

#include "qsc.hpp"

namespace qsc {

  typedef unsigned long long int big;
  
  const std::string SCAN_OPTION_LINEAR = "linear";
  const std::string SCAN_OPTION_LOG = "log";
  const std::string SCAN_OPTION_2_SIDED_LOG = "2 sided log";
  const std::string SCAN_OPTION_2_SIDED_LOG_EXCEPT_Z0s1 = "2 sided log except Z0s1";
  
  class Scan {
  private:
    void defaults();
    
  public:
    Qsc q;
    std::string eta_bar_scan_option, sigma0_scan_option;
    std::string B2s_scan_option, B2c_scan_option, fourier_scan_option;
    qscfloat eta_bar_min, eta_bar_max, sigma0_min, sigma0_max;
    qscfloat B2s_min, B2s_max, B2c_min, B2c_max;
    Vector R0c_min, R0c_max, R0s_min, R0s_max, Z0c_min, Z0c_max, Z0s_min, Z0s_max;
    qscfloat max_seconds;
    big n_scan, max_keep_per_proc;
    qscfloat min_iota_to_keep, min_R0_to_keep, max_elongation_to_keep;
    qscfloat min_L_grad_B_to_keep, min_L_grad_grad_B_to_keep, min_r_singularity_to_keep;
    qscfloat max_d2_volume_d_psi2_to_keep, min_DMerc_to_keep;
    bool keep_all;

    Scan();
    void run(std::string);
    void input(std::string);
    
  };
}

#endif

