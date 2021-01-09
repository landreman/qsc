#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include "opt.hpp"

using namespace qsc;

void Opt::optimize() {
  std::cout << "optimizing..." << std::endl;
  n_terms = 4;
  n_parameters = 4;
  
  gsl_vector *gsl_residual = gsl_vector_alloc(n_terms);
  gsl_vector *gsl_state_vector = gsl_vector_alloc(n_parameters);
  gsl_multifit_nlinear_fdf gsl_optimizer;
  gsl_multifit_nlinear_parameters gsl_optimizer_params = gsl_multifit_nlinear_default_parameters();
}
