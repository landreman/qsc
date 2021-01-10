#include <string>
#include <stdexcept>
#include "qsc.hpp"

using namespace qsc;

/** Initialize the parameters for a Qsc object from one of the
 *  configurations discussed in our published papers.
 */
Qsc::Qsc(std::string config_name) {
  defaults();
  
  if (config_name.compare("r1 section 5.1") == 0) {
    nfp = 3;
    eta_bar = -0.9;
    order_r_option = ORDER_R_OPTION_R1;

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
    order_r_option = ORDER_R_OPTION_R1;

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
    order_r_option = ORDER_R_OPTION_R1;

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
    eta_bar = 0.64;
    B2c = -0.00322;
    order_r_option = ORDER_R_OPTION_R2p1;

    R0c.resize(3, 0.0);
    R0s.resize(3, 0.0);
    Z0c.resize(3, 0.0);
    Z0s.resize(3, 0.0);
    R0c[0] = 1.0;
    R0c[1] = 0.155;
    R0c[2] = 0.0102;
    Z0s[1] = 0.154;
    Z0s[2] = 0.0111;

  } else if (config_name.compare("r2 section 5.2") == 0 ||
	     config_name.compare("2") == 0) {
    nfp = 2;
    eta_bar = 0.632;
    B2c = -0.158;
    order_r_option = ORDER_R_OPTION_R2p1;

    R0c.resize(4, 0.0);
    R0s.resize(4, 0.0);
    Z0c.resize(4, 0.0);
    Z0s.resize(4, 0.0);
    R0c[0] = 1.0;
    R0c[1] = 0.173;
    R0c[2] = 0.0168;
    R0c[3] = 0.00101;
    
    Z0s[1] = 0.159;
    Z0s[2] = 0.0165;
    Z0s[3] = 0.000985;

  } else if (config_name.compare("r2 section 5.3") == 0 ||
	     config_name.compare("3") == 0) {
    nfp = 2;
    eta_bar = 0.95;
    B2c = -0.7;
    I2 = 0.9;
    p2 = -600000.0;
    order_r_option = ORDER_R_OPTION_R2p1;
    
    R0c.resize(2, 0.0);
    R0s.resize(2, 0.0);
    Z0c.resize(2, 0.0);
    Z0s.resize(2, 0.0);
    R0c[0] = 1.0;
    R0c[1] = 0.09;
    Z0s[1] = -0.09;

  } else if (config_name.compare("r2 section 5.4") == 0 ||
	     config_name.compare("4") == 0) {
    nfp = 4;
    eta_bar = 1.569;
    B2c = 0.1348;
    order_r_option = ORDER_R_OPTION_R2p1;

    R0c.resize(5, 0.0);
    R0s.resize(5, 0.0);
    Z0c.resize(5, 0.0);
    Z0s.resize(5, 0.0);
    R0c[0] = 1.0;
    R0c[1] = 0.17;
    R0c[2] = 0.01804;
    R0c[3] = 0.001409;
    R0c[4] = 5.877e-05;
    
    Z0s[1] = 0.1581;
    Z0s[2] = 0.01820;
    Z0s[3] = 0.001548;
    Z0s[4] = 7.772e-05;

  } else if (config_name.compare("r2 section 5.5") == 0 ||
	     config_name.compare("5") == 0) {
    nfp = 5;
    eta_bar = 2.5;
    sigma0 = 0.3;
    B2c = 1.0;
    B2s = 3.0;
    I2 = 1.6;
    p2 = -0.5e7;
    order_r_option = ORDER_R_OPTION_R2p1;

    R0c.resize(2, 0.0);
    R0s.resize(2, 0.0);
    Z0c.resize(2, 0.0);
    Z0s.resize(2, 0.0);
    R0c[0] = 1.0;
    R0c[1] = 0.3;
    Z0s[1] = 0.3;

  } else {
    throw std::runtime_error("Configuration name not recognized");
  }

  init();
  calculate();
}
