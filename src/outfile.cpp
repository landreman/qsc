#include <iostream>
#include <stdexcept>
#include "qsc.hpp"

std::string qsc::outfile(std::string directory_and_infile) {
  size_t slash_pos = directory_and_infile.rfind('/'); // Find the last appearance of '/'
  std::string directory = "";
  std::string infile(directory_and_infile);
  if (slash_pos != std::string::npos) {
    // We did find '/'
    directory = directory_and_infile.substr(0, slash_pos + 1);
    infile = directory_and_infile.substr(slash_pos + 1, directory_and_infile.size() - slash_pos - 1);
  }
  
  if (infile.size() < 7) {
    throw std::runtime_error("Input file must begin with 'qsc_in.'");
  }
  std::string infile_start = infile.substr(0, 7);
  if (infile_start != "qsc_in.") {
    throw std::runtime_error("Input file must begin with 'qsc_in.'");
  }
  std::string directory_and_outfile = directory + "qsc_out." + infile.substr(7, infile.size() - 7) + ".nc";
  return directory_and_outfile;
}
