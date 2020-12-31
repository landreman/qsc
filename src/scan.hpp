#ifndef QSC_SCAN_H
#define QSC_SCAN_H

#include <mpi.h>
#include "qsc.hpp"

#ifdef SINGLE
#define MPI_QSCFLOAT MPI_FLOAT
#else
#define MPI_QSCFLOAT MPI_DOUBLE
#endif

namespace qsc {

  typedef unsigned long int big;
  
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
    big n_scan, attempts;
    int max_keep_per_proc, max_attempts_per_proc; // Can I read in a "big" from toml?
    qscfloat min_R0_to_keep, min_iota_to_keep, max_elongation_to_keep;
    qscfloat min_L_grad_B_to_keep, min_L_grad_grad_B_to_keep;
    qscfloat max_B20_variation_to_keep, min_r_singularity_to_keep;
    qscfloat max_d2_volume_d_psi2_to_keep, min_DMerc_to_keep;
    bool keep_all, deterministic;
    big rejected_due_to_R0_crude, rejected_due_to_R0, rejected_due_to_curvature;
    big rejected_due_to_iota, rejected_due_to_elongation;
    big rejected_due_to_L_grad_B, rejected_due_to_L_grad_grad_B;
    big rejected_due_to_B20_variation, rejected_due_to_r_singularity;
    big rejected_due_to_d2_volume_d_psi2, rejected_due_to_DMerc;

    Vector scan_eta_bar, scan_sigma0, scan_B2s, scan_B2c;
    Matrix scan_R0c, scan_R0s, scan_Z0c, scan_Z0s;
    Vector scan_min_R0, scan_max_curvature;
    Vector scan_iota, scan_max_elongation;
    Vector scan_min_L_grad_B, scan_min_L_grad_grad_B;
    Vector scan_r_singularity, scan_B20_variation;
    Vector scan_d2_volume_d_psi2, scan_DMerc_times_r2;
    std::valarray<int> scan_helicity;
    
    Scan();
    void run(std::string);
    void input(std::string);
    void random();
  };
}

#endif

