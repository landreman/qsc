#ifndef QSC_OPT_H
#define QSC_OPT_H

#include <valarray>
#include "qsc.hpp"

namespace qsc {

  class Opt {
  private:
    void defaults();
    
  public:
    Qsc q;
    int max_iter, n_iter;
    int n_parameters, n_terms;
    int verbose;
    std::string outfilename;
    Vector sqrt_d_l_d_phi;

    bool vary_eta_bar, vary_sigma0;
    bool vary_B2c, vary_B2s;
    std::valarray<bool> vary_R0c, vary_R0s, vary_Z0c, vary_Z0s;
    
    qscfloat weight_B20, weight_iota, target_iota;
    qscfloat weight_R0, min_R0;
    qscfloat weight_d2_volume_d_psi2, max_d2_volume_d_psi2;
    qscfloat weight_XY2, weight_XY2Prime;
    qscfloat weight_XY3, weight_XY3Prime;
    qscfloat weight_grad_grad_B;

    Vector iter_objective_function;
    Vector iter_B20_term, iter_iota_term;
    Vector iter_R0_term, iter_d2_volume_d_psi2_term;
    Vector iter_XY2_term, iter_XY2Prime_term;
    Vector iter_XY3_term, iter_XY3Prime_term;
    Vector iter_grad_grad_B_term;
    
    Vector iter_eta_bar, iter_sigma0, iter_B2s, iter_B2c;
    Matrix iter_R0c, iter_R0s, iter_Z0c, iter_Z0s;
    Vector iter_min_R0, iter_max_curvature;
    Vector iter_iota, iter_max_elongation;
    Vector iter_min_L_grad_B, iter_min_L_grad_grad_B;
    Vector iter_r_singularity, iter_B20_variation, iter_B20_residual;
    Vector iter_d2_volume_d_psi2, iter_DMerc_times_r2;
    Vector iter_standard_deviation_of_R, iter_standard_deviation_of_Z;
    
    Opt();
    void run(std::string);
    void allocate();
    void input(std::string);
    void optimize();
    void write_netcdf();
  };
}

#endif

