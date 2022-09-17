#include "opt.hpp"

using namespace qsc;

void Opt::defaults() {
  // Set defaults.
  verbose = 1;
  max_iter = 3000;
  make_names = true;
  algorithm = GSL_LM;
  fourier_refine = 0;
  refine_angle_shift = false;
  toml_group = "opt";
  diff_method = DIFF_METHOD_FORWARD;
  n_evals = 0;

  vary_eta_bar = true;
  vary_sigma0 = false;
  vary_B2c = true;
  vary_B2s = false;

  // These defaults are also set in Opt::input()
  vary_R0c.resize(q.R0c.size(), true);
  vary_R0s.resize(q.R0c.size(), false);
  vary_Z0c.resize(q.R0c.size(), false);
  vary_Z0s.resize(q.R0c.size(), true);
  vary_fc.resize(q.R0c.size(), false);
  vary_fs.resize(q.R0c.size(), false);
  vary_R0c[0] = false;
  vary_Z0s[0] = false;
  
  weight_B20 = -1.0;
  weight_iota = -1.0;
  target_iota = 0.4;
  weight_elongation = -1.0;
  weight_curvature = -1.0;
  weight_R0 = -1.0;
  min_R0 = 0.3;
  weight_d2_volume_d_psi2 = -1.0;
  max_d2_volume_d_psi2 = 0.0;
  weight_DMerc_times_r2 = -1.0;
  min_DMerc_times_r2 = 0.0;
  weight_XY2 = -1.0;
  weight_XY2Prime = -1.0;
  weight_XY2PrimePrime = -1.0;
  weight_Z2 = -1.0;
  weight_Z2Prime = -1.0;
  weight_XY3 = -1.0;
  weight_XY3Prime = -1.0;
  weight_XY3PrimePrime = -1.0;
  weight_grad_B = -1.0;
  weight_grad_grad_B = -1.0;
  weight_r_singularity = -1.0;
  weight_axis_length = -1.0;
  target_axis_length = 0.0;
  weight_standard_deviation_of_R = -1.0;
  weight_B20_mean = -1.0;
}

Opt::Opt() {
  defaults();
}

void Opt::run(std::string directory_and_infile) {
  outfilename = qsc::outfile(directory_and_infile);
  input(directory_and_infile);
  allocate();
  optimize();
  write_netcdf();
}
