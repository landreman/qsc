#include <string>
#include <stdexcept>
#include "qsc.hpp"

using namespace qsc;

/** Initialize the parameters for a Qsc object from one of the
 *  configurations discussed in our published papers.
 */
Qsc::Qsc(std::string config_name) :
  // Call constructor of member objects:
  d_d_phi(1, 1),
  d_d_varphi(1, 1)
{
  // Defaults common to most or all of these configurations:
  verbose = 1;
  nphi = 31;
  sG = 1;
  spsi = 1;
  B0 = 1.0;
  I2 = 0.0;
  sigma0 = 0.0;
  B2s = 0.0;
  B2c = 0.0;
  
  if (config_name.compare("r1 section 5.1") == 0) {
    nfp = 3;
    eta_bar = -0.9;

    R0c.resize(2, 0.0);
    R0s.resize(2, 0.0);
    Z0c.resize(2, 0.0);
    Z0s.resize(2, 0.0);
    R0c[0] = 1.0;
    R0c[1] = 0.045;
    Z0s[1] = -0.045;
    
  } else if (config_name.compare("r1 section 5.2") == 0) {
    nfp = 4;
    eta_bar = -2.25;

    R0c.resize(2, 0.0);
    R0s.resize(2, 0.0);
    Z0c.resize(2, 0.0);
    Z0s.resize(2, 0.0);
    R0c[0] = 1.0;
    R0c[1] = 0.265;
    Z0s[1] = -0.21;

  } else if (config_name.compare("r1 section 5.3") == 0) {
    nfp = 3;
    eta_bar = -1.1;
    sigma0 = -0.6;

    R0c.resize(2, 0.0);
    R0s.resize(2, 0.0);
    Z0c.resize(2, 0.0);
    Z0s.resize(2, 0.0);
    R0c[0] = 1.0;
    R0c[1] = 0.042;
    Z0s[1] = -0.042;
    Z0c[1] = -0.025;

  } else if (config_name.compare("r2 section 5.1") == 0 ||
	     config_name.compare("1") == 0) {
    nfp = 2;
    eta_bar = -0.64;
    B2c = -0.00322;

    R0c.resize(3, 0.0);
    R0s.resize(3, 0.0);
    Z0c.resize(3, 0.0);
    Z0s.resize(3, 0.0);
    R0c[0] = 1.0;
    R0c[1] = 0.155;
    R0c[2] = 0.0102;
    Z0s[1] = 0.154;
    Z0s[2] = 0.0111;

  } else {
    throw std::runtime_error("Configuration name not recognized");
  }

  calculate();
}
