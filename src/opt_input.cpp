#include <stdexcept>
#include "qsc.hpp"
#include "opt.hpp"
#include "toml.hpp"
#include "toml_util.hpp"

using namespace qsc;

/** Read in an optimization input file
 */
void Opt::input(std::string filename) {
  int j;

  q.input(filename);
  // Now that we know the size of the Qsc axis shape arrays,
  // set defaults for the vary_* arrays:
  vary_R0c.resize(q.R0c.size(), true);
  vary_R0s.resize(q.R0c.size(), false);
  vary_Z0c.resize(q.R0c.size(), false);
  vary_Z0s.resize(q.R0c.size(), true);
  vary_R0c[0] = false;
  vary_Z0s[0] = false;
    
  auto toml_file = toml::parse(filename);
  auto indata = toml::find(toml_file, toml_group);
  
  std::vector<std::string> varlist;

  toml_read(varlist, indata, "vary_eta_bar", vary_eta_bar);
  toml_read(varlist, indata, "vary_sigma0", vary_sigma0);
  toml_read(varlist, indata, "vary_B2c", vary_B2c);
  toml_read(varlist, indata, "vary_B2s", vary_B2s);
  toml_read(varlist, indata, "vary_R0c", vary_R0c);
  toml_read(varlist, indata, "vary_R0s", vary_R0s);
  toml_read(varlist, indata, "vary_Z0c", vary_Z0c);
  toml_read(varlist, indata, "vary_Z0s", vary_Z0s);
  
  toml_read(varlist, indata, "weight_B20", weight_B20);
  toml_read(varlist, indata, "weight_iota", weight_iota);
  toml_read(varlist, indata, "target_iota", target_iota);
  toml_read(varlist, indata, "weight_elongation", weight_elongation);
  toml_read(varlist, indata, "weight_curvature", weight_curvature);
  toml_read(varlist, indata, "weight_R0", weight_R0);
  toml_read(varlist, indata, "min_R0", min_R0);
  toml_read(varlist, indata, "weight_d2_volume_d_psi2", weight_d2_volume_d_psi2);
  toml_read(varlist, indata, "max_d2_volume_d_psi2", max_d2_volume_d_psi2);
  toml_read(varlist, indata, "weight_XY2", weight_XY2);
  toml_read(varlist, indata, "weight_XY2Prime", weight_XY2Prime);
  toml_read(varlist, indata, "weight_Z2", weight_Z2);
  toml_read(varlist, indata, "weight_Z2Prime", weight_Z2Prime);
  toml_read(varlist, indata, "weight_XY3", weight_XY3);
  toml_read(varlist, indata, "weight_XY3Prime", weight_XY3Prime);
  toml_read(varlist, indata, "weight_grad_B", weight_grad_B);
  toml_read(varlist, indata, "weight_grad_grad_B", weight_grad_grad_B);
  toml_read(varlist, indata, "weight_r_singularity", weight_r_singularity);
  toml_read(varlist, indata, "weight_axis_length", weight_axis_length);
  toml_read(varlist, indata, "weight_standard_deviation_of_R", weight_standard_deviation_of_R);

  toml_read(varlist, indata, "max_iter", max_iter);
  toml_read(varlist, indata, "verbose", verbose);
  toml_read(varlist, indata, "make_names", make_names);
  toml_read(varlist, indata, "fourier_refine", fourier_refine);

  std::string algorithm_str = "";
  toml_read(varlist, indata, "algorithm", algorithm_str);
  if (algorithm_str.compare("lm") == 0) {
    algorithm = GSL_LM;
  } else if (algorithm_str.compare("dogleg") == 0) {
    algorithm = GSL_DOGLEG;
  } else if (algorithm_str.compare("ddogleg") == 0) {
    algorithm = GSL_DDOGLEG;
  } else if (algorithm_str.compare("subspace2d") == 0) {
    algorithm = GSL_SUBSPACE2D;
  } else if (algorithm_str.compare("") == 0) {
    // Stick with the default
  } else {
    throw std::runtime_error("Unrecognized algorithm");
  }

  std::string diff_method_str = "";
  toml_read(varlist, indata, "diff_method", diff_method_str);
  if (diff_method_str.compare("forward") == 0) {
    diff_method = DIFF_METHOD_FORWARD;
  } else if (diff_method_str.compare("centered") == 0) {
    diff_method = DIFF_METHOD_CENTERED;
  } else if (diff_method_str.compare("") == 0) {
    // Stick with the default
  } else {
    throw std::runtime_error("Unrecognized diff_method");
  }
  
  toml_unused(varlist, indata);
  
  if (verbose > 0) {
    std::cout << "----- Optimization parameters -----" << std::endl;
    std::cout << "max_iter: " << max_iter << std::endl;
    std::cout << "fourier_refine: " << fourier_refine << std::endl;
    std::cout << "algorithm: ";
    switch (algorithm) {
    case GSL_LM:
      std::cout << "Levenberg-Marquardt";
      break;
    case GSL_DOGLEG:
      std::cout << "Dogleg";
      break;
    case GSL_DDOGLEG:
      std::cout << "Double dogleg";
      break;
    case GSL_SUBSPACE2D:
      std::cout << "Subspace 2D";
      break;
    default:
      throw std::runtime_error("Should not get here");
    }
    std::cout << std::endl;
    
    std::cout << "vary_eta_bar: " << vary_eta_bar << std::endl;
    std::cout << "vary_sigma0: " << vary_sigma0 << std::endl;
    std::cout << "vary_B2c: " << vary_B2c << std::endl;
    std::cout << "vary_B2s: " << vary_B2s << std::endl;
    
    std::cout << "vary_R0c:";
    for (j = 0; j < vary_R0c.size(); j++) std::cout << " " << vary_R0c[j];
    std::cout << std::endl;
    
    std::cout << "vary_R0s:";
    for (j = 0; j < vary_R0s.size(); j++) std::cout << " " << vary_R0s[j];
    std::cout << std::endl;
    
    std::cout << "vary_Z0c:";
    for (j = 0; j < vary_Z0c.size(); j++) std::cout << " " << vary_Z0c[j];
    std::cout << std::endl;
    
    std::cout << "vary_Z0s:";
    for (j = 0; j < vary_Z0s.size(); j++) std::cout << " " << vary_Z0s[j];
    std::cout << std::endl;
    
    std::cout << "weight_B20: " << weight_B20 << std::endl;
    std::cout << "weight_iota: " << weight_iota << std::endl;
    std::cout << "target_iota: " << target_iota << std::endl;
    std::cout << "weight_elongation: " << weight_elongation << std::endl;
    std::cout << "weight_curvature: " << weight_curvature << std::endl;
    std::cout << "weight_R0: " << weight_R0 << std::endl;
    std::cout << "min_R0: " << min_R0 << std::endl;
    std::cout << "weight_d2_volume_d_psi2: " << weight_d2_volume_d_psi2 << std::endl;
    std::cout << "max_d2_volume_d_psi2: " << max_d2_volume_d_psi2 << std::endl;
    std::cout << "weight_XY2: " << weight_XY2 << std::endl;
    std::cout << "weight_XY2Prime: " << weight_XY2Prime << std::endl;
    std::cout << "weight_Z2: " << weight_Z2 << std::endl;
    std::cout << "weight_Z2Prime: " << weight_Z2Prime << std::endl;
    std::cout << "weight_XY3: " << weight_XY3 << std::endl;
    std::cout << "weight_XY3Prime: " << weight_XY3Prime << std::endl;
    std::cout << "weight_grad_B: " << weight_grad_B << std::endl;
    std::cout << "weight_grad_grad_B: " << weight_grad_grad_B << std::endl;
    std::cout << "weight_r_singularity: " << weight_r_singularity << std::endl;
  }
}
