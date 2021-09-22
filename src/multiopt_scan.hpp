#ifndef QSC_MULTIOPT_SCAN_H
#define QSC_MULTIOPT_SCAN_H

#include <vector>
#include <valarray>
#include <mpi.h>
#include "multiopt.hpp"

namespace qsc {

  class MultiOptScan {
  private:
    MultiOpt mo, mo_ref;
    void defaults();
    
    enum {ATTEMPTS,
      KEPT,
      REJECTED_DUE_TO_R0,
      REJECTED_DUE_TO_CURVATURE,
      REJECTED_DUE_TO_IOTA,
      REJECTED_DUE_TO_ELONGATION,
      REJECTED_DUE_TO_L_GRAD_B,
      REJECTED_DUE_TO_B20_VARIATION,
      REJECTED_DUE_TO_L_GRAD_GRAD_B,
      REJECTED_DUE_TO_D2_VOLUME_D_PSI2,
      REJECTED_DUE_TO_DMERC,
      REJECTED_DUE_TO_R_SINGULARITY,
      N_FILTERS};

  public:
    MPI_Comm mpi_comm;
    qscfloat max_seconds, save_period;
    big n_scan, n_scan_all, filters[N_FILTERS];
    qscfloat filter_fractions[N_FILTERS];
    int max_keep_per_proc, max_attempts_per_proc; // Can I read in a "big" from toml?
    qscfloat min_R0_to_keep, min_iota_to_keep, max_elongation_to_keep;
    qscfloat min_L_grad_B_to_keep, min_L_grad_grad_B_to_keep;
    qscfloat max_B20_variation_to_keep, min_r_singularity_to_keep;
    qscfloat max_d2_volume_d_psi2_to_keep, min_DMerc_to_keep;
    bool keep_all;
    int verbose;
    std::string outfilename;
    int ndim;
    std::vector<std::string> params;
    Vector params_max, params_min;
    std::valarray<bool> params_log;
    std::valarray<int> params_n;
    std::vector<Vector> params_vals;
    
    MultiOptScan();
    void run(std::string);
    void input(std::string);
    void init();
    void scan();
    //void write_netcdf();
  };
}

#endif

