#include <string>
#include <stdexcept>
#include "multiopt.hpp"
#include "toml.hpp"
#include "toml_util.hpp"

using namespace qsc;

void MultiOpt::defaults() {
  // Set defaults.
  verbose = 0;
}

MultiOpt::MultiOpt() {
  defaults();
}

void MultiOpt::run(std::string directory_and_infile) {
  std::string outfile_base_name = qsc::outfile(directory_and_infile);
  input(directory_and_infile);
  for (int jopt = 0; jopt < opts.size(); jopt++) {
    opts[jopt].outfilename = outfile_base_name.substr(0, outfile_base_name.size() - 3)
      + "_" + std::to_string(jopt) + ".nc";
    
    if (verbose > 0) std::cout << "Output filename " << jopt << " will be " << opts[jopt].outfilename << std::endl;
  }
  
  optimize();
  write_netcdf();
}

/**
 * Read a QSC input file for a MultiOpt job.
 */
void MultiOpt::input(std::string filename) {
  auto toml_file = toml::parse(filename);
  auto indata = toml::find(toml_file, "multiopt");

  std::vector<std::string> varlist;
  int nopts;
  toml_read(varlist, indata, "nopts", nopts);
  toml_read(varlist, indata, "verbose", verbose);
  toml_unused(varlist, indata);

  if (nopts < 1) throw std::runtime_error("In the multiopt parameter list, nopts must be at least 1.");
  opts.resize(nopts);
  for (int jopt = 0; jopt < nopts; jopt++) {
    opts[jopt].make_names = false; // Do not save a file with the names of the residuals.
    opts[jopt].toml_group = "opt" + std::to_string(jopt);
    if (verbose > 0) std::cout << "About to read toml group " << opts[jopt].toml_group << std::endl;
    opts[jopt].input(filename);
  }
}

/**
 * Carry out the multi-stage optimization.
 */
void MultiOpt::optimize() {
  for (int jopt = 0; jopt < opts.size(); jopt++) {
    if (jopt > 0) {
      // Transfer data from the end of the previous optimization stage to the next one.
      // opts[jopt].q = Qsc();

      if (verbose > 0) std::cout << "Transferring data to optimization stage " << jopt << " from previous stage." << std::endl;
      /*
      opts[jopt].q.nfp = opts[jopt - 1].q.nfp;
      opts[jopt].q.nphi = opts[jopt - 1].q.nphi;
      opts[jopt].q.order_r_option = opts[jopt - 1].q.order_r_option;
      opts[jopt].q.eta_bar = opts[jopt - 1].q.eta_bar;
      opts[jopt].q.sigma0 = opts[jopt - 1].q.sigma0;
      opts[jopt].q.B2c = opts[jopt - 1].q.B2c;
      opts[jopt].q.B2s = opts[jopt - 1].q.B2s;
      if (verbose > 0) std::cout << "Size of R0c before transfer: " << opts[jopt].q.R0c.size() << std::endl;
      opts[jopt].q.R0c = opts[jopt - 1].q.R0c;
      if (verbose > 0) std::cout << "Size of R0c after transfer: " << opts[jopt].q.R0c.size() << std::endl;
      opts[jopt].q.R0s = opts[jopt - 1].q.R0s;
      opts[jopt].q.Z0c = opts[jopt - 1].q.Z0c;
      opts[jopt].q.Z0s = opts[jopt - 1].q.Z0s;
      opts[jopt].q.I2 = opts[jopt - 1].q.I2;
      opts[jopt].q.p2 = opts[jopt - 1].q.p2;
      */
      opts[jopt].q = opts[jopt - 1].q;
    }

    // Run the given stage of the optimization.
    if (verbose > 0) {
      std::cout << "################################################" << std::endl;
      std::cout << "######## Beginning optimization stage " << jopt << " ########" << std::endl;
      std::cout << "################################################" << std::endl;
    }
    opts[jopt].allocate();
    opts[jopt].optimize();
    if (verbose > 0) {
      std::cout << "################################################" << std::endl;
      std::cout << "######## Done with optimization stage " << jopt << " ########" << std::endl;
      std::cout << "################################################" << std::endl;
    }
  }
}

void MultiOpt::write_netcdf() {
  for (int jopt = 0; jopt < opts.size(); jopt++) {
    opts[jopt].write_netcdf();
  }
}
