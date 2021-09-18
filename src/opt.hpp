#ifndef QSC_OPT_H
#define QSC_OPT_H

#include <valarray>
#include <vector>
#include <gsl/gsl_vector.h>
#include "qsc.hpp"

namespace qsc {

  typedef enum {
    GSL_LM,
    GSL_DOGLEG,
    GSL_DDOGLEG,
    GSL_SUBSPACE2D
  } algorithm_type;

  typedef enum {
    DIFF_METHOD_FORWARD,
    DIFF_METHOD_CENTERED
  } diff_method_type;
  
  class Opt {
  private:    
    void defaults();
    
  public:
    Qsc q;
    int max_iter, n_iter;
    int j_fourier_refine;
    int n_parameters, n_terms;
    int verbose;
    std::string outfilename;
    std::string toml_group;
    Vector arclength_factor;
    bool make_names;
    std::vector<std::string> state_vector_names, residual_names;
    algorithm_type algorithm;
    Vector residuals;
    int fourier_refine;
    diff_method_type diff_method;

    bool vary_eta_bar, vary_sigma0;
    bool vary_B2c, vary_B2s;
    std::valarray<bool> vary_R0c, vary_R0s, vary_Z0c, vary_Z0s;
    
    qscfloat weight_B20, weight_iota, target_iota;
    qscfloat weight_elongation, weight_curvature;
    qscfloat weight_R0, min_R0;
    qscfloat weight_d2_volume_d_psi2, max_d2_volume_d_psi2;
    qscfloat weight_XY2, weight_XY2Prime;
    qscfloat weight_Z2, weight_Z2Prime;
    qscfloat weight_XY3, weight_XY3Prime;
    qscfloat weight_grad_B, weight_grad_grad_B, weight_r_singularity;

    qscfloat objective_function;
    qscfloat B20_term, iota_term;
    qscfloat elongation_term, curvature_term;
    qscfloat R0_term, d2_volume_d_psi2_term;
    qscfloat XY2_term, XY2Prime_term;
    qscfloat Z2_term, Z2Prime_term;
    qscfloat XY3_term, XY3Prime_term;
    qscfloat grad_B_term, grad_grad_B_term, r_singularity_term;

    Vector iter_objective_function;
    Vector iter_B20_term, iter_iota_term;
    Vector iter_elongation_term, iter_curvature_term;
    Vector iter_R0_term, iter_d2_volume_d_psi2_term;
    Vector iter_XY2_term, iter_XY2Prime_term;
    Vector iter_Z2_term, iter_Z2Prime_term;
    Vector iter_XY3_term, iter_XY3Prime_term;
    Vector iter_grad_B_term, iter_grad_grad_B_term, iter_r_singularity_term;
    
    Vector iter_eta_bar, iter_sigma0, iter_B2s, iter_B2c;
    Matrix iter_R0c, iter_R0s, iter_Z0c, iter_Z0s;
    Vector iter_min_R0, iter_max_curvature;
    Vector iter_iota, iter_max_elongation;
    Vector iter_min_L_grad_B, iter_min_L_grad_grad_B;
    Vector iter_r_singularity, iter_B20_variation, iter_B20_residual;
    Vector iter_d2_volume_d_psi2, iter_DMerc_times_r2;
    Vector iter_standard_deviation_of_R, iter_standard_deviation_of_Z;
    std::valarray<int> iter_fourier_refine_step;
    
    Opt();
    void run(std::string);
    void allocate();
    void init_parameters();
    void init_residuals();
    void input(std::string);
    void optimize();
    void write_netcdf();
    void set_state_vector(qscfloat*);
    void unpack_state_vector(qscfloat*);
    void set_residuals(gsl_vector *);
  };
}

#endif

