#include "qsc.hpp"
#include "scan.hpp"
#include "toml.hpp"
#include "toml_util.hpp"

using namespace qsc;

/** Read in a scan input file
 */
void Scan::input(std::string filename) {

  q.input(filename);
    
  auto toml_file = toml::parse(filename);
  auto indata = toml::find(toml_file, "scan");
  
  std::vector<std::string> varlist;
  
  toml_read(varlist, indata, "eta_bar_scan_option", eta_bar_scan_option);
  toml_read(varlist, indata, "sigma0_scan_option", sigma0_scan_option);
  toml_read(varlist, indata, "B2c_scan_option", B2c_scan_option);
  toml_read(varlist, indata, "B2s_scan_option", B2s_scan_option);
  toml_read(varlist, indata, "fourier_scan_option", fourier_scan_option);
  
  toml_read(varlist, indata, "eta_bar_min", eta_bar_min);
  toml_read(varlist, indata, "eta_bar_max", eta_bar_max);
  toml_read(varlist, indata, "sigma0_min", sigma0_min);
  toml_read(varlist, indata, "sigma0_max", sigma0_max);
  toml_read(varlist, indata, "B2c_min", B2c_min);
  toml_read(varlist, indata, "B2c_max", B2c_max);
  toml_read(varlist, indata, "B2s_min", B2s_min);
  toml_read(varlist, indata, "B2s_max", B2s_max);

  toml_read(varlist, indata, "R0c_min", R0c_min);
  toml_read(varlist, indata, "R0c_max", R0c_max);
  toml_read(varlist, indata, "R0s_min", R0s_min);
  toml_read(varlist, indata, "R0s_max", R0s_max);
  toml_read(varlist, indata, "Z0c_min", Z0c_min);
  toml_read(varlist, indata, "Z0c_max", Z0c_max);
  toml_read(varlist, indata, "Z0s_min", Z0s_min);
  toml_read(varlist, indata, "Z0s_max", Z0s_max);

  toml_read(varlist, indata, "deterministic", deterministic);
  toml_read(varlist, indata, "max_seconds", max_seconds);
  toml_read(varlist, indata, "max_keep_per_proc", max_keep_per_proc);
  toml_read(varlist, indata, "max_attempts_per_proc", max_attempts_per_proc);

  toml_read(varlist, indata, "keep_all", keep_all);
  toml_read(varlist, indata, "min_R0_to_keep", min_R0_to_keep);
  toml_read(varlist, indata, "min_iota_to_keep", min_iota_to_keep);
  toml_read(varlist, indata, "max_elongation_to_keep", max_elongation_to_keep);
  toml_read(varlist, indata, "min_L_grad_B_to_keep", min_L_grad_B_to_keep);
  toml_read(varlist, indata, "min_L_grad_grad_B_to_keep", min_L_grad_grad_B_to_keep);
  toml_read(varlist, indata, "max_B20_variation_to_keep", max_B20_variation_to_keep);
  toml_read(varlist, indata, "min_r_singularity_to_keep", min_r_singularity_to_keep);
  toml_read(varlist, indata, "min_DMerc_to_keep", min_DMerc_to_keep);
  toml_read(varlist, indata, "max_d2_volume_d_psi2_to_keep", max_d2_volume_d_psi2_to_keep);

  toml_unused(varlist, indata);
  
  // Now, make {R0c, R0s, Z0c, Z0s} have the same size.
  std::size_t newsize = 0;
  newsize = std::max(newsize, R0c_min.size());
  newsize = std::max(newsize, R0s_min.size());
  newsize = std::max(newsize, Z0c_min.size());
  newsize = std::max(newsize, Z0s_min.size());
  newsize = std::max(newsize, R0c_max.size());
  newsize = std::max(newsize, R0s_max.size());
  newsize = std::max(newsize, Z0c_max.size());
  newsize = std::max(newsize, Z0s_max.size());
  pad_vector(R0c_min, newsize);
  pad_vector(R0s_min, newsize);
  pad_vector(Z0c_min, newsize);
  pad_vector(Z0s_min, newsize);
  pad_vector(R0c_max, newsize);
  pad_vector(R0s_max, newsize);
  pad_vector(Z0c_max, newsize);
  pad_vector(Z0s_max, newsize);

  std::cout << "----- Scan parameters -----" << std::endl;
  std::cout << "eta_bar range: " << eta_bar_min << " to " << eta_bar_max << std::endl;
  std::cout << "sigma0 range: " << sigma0_min << " to " << sigma0_max << std::endl;
  std::cout << "B2c range: " << B2c_min << " to " << B2c_max << std::endl;
  std::cout << "B2s range: " << B2s_min << " to " << B2s_max << std::endl;

  std::cout << "R0c_min: " << R0c_min << std::endl;
  std::cout << "R0c_max: " << R0c_max << std::endl;
  std::cout << "R0s_min: " << R0s_min << std::endl;
  std::cout << "R0s_max: " << R0s_max << std::endl;
  std::cout << "Z0c_min: " << Z0c_min << std::endl;
  std::cout << "Z0c_max: " << Z0c_max << std::endl;
  std::cout << "Z0s_min: " << Z0s_min << std::endl;
  std::cout << "Z0s_max: " << Z0s_max << std::endl;
  
  std::cout << "eta_bar_scan_option: " << eta_bar_scan_option << std::endl;
  std::cout << "sigma0_scan_option: " << sigma0_scan_option << std::endl;
  std::cout << "B2c_scan_option: " << B2c_scan_option << std::endl;
  std::cout << "B2s_scan_option: " << B2s_scan_option << std::endl;
  std::cout << "fourier_scan_option: " << fourier_scan_option << std::endl;

  std::cout << "max_seconds: " << max_seconds << std::endl;
  std::cout << "max_keep_per_proc: " << max_keep_per_proc << std::endl;
  std::cout << "deterministic: " << deterministic << std::endl;
  std::cout << "keep_all: " << keep_all << std::endl;
  if (!keep_all) {
    std::cout << "min_R0_to_keep: " << min_R0_to_keep << std::endl;
    std::cout << "min_iota_to_keep: " << min_iota_to_keep << std::endl;
    std::cout << "max_elongation_to_keep: " << max_elongation_to_keep << std::endl;
    std::cout << "min_L_grad_B_to_keep: " << min_L_grad_B_to_keep << std::endl;
    std::cout << "min_L_grad_grad_B_to_keep: " << min_L_grad_grad_B_to_keep << std::endl;
    std::cout << "max_B20_variation_to_keep: " << max_B20_variation_to_keep << std::endl;
    std::cout << "min_r_singularity_to_keep: " << min_r_singularity_to_keep << std::endl;
    std::cout << "max_d2_volume_d_psi2_to_keep: " << max_d2_volume_d_psi2_to_keep << std::endl;
    std::cout << "min_DMerc_to_keep: " << min_DMerc_to_keep << std::endl;
  }
}
