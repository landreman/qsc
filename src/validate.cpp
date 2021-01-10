#include <string>
#include <stdexcept>
#include "qsc.hpp"

using namespace qsc;

/** Handle string input options.
 */
void Qsc::validate() {
  
  if (order_r_option.compare(ORDER_R_OPTION_R1) == 0) {
    at_least_order_r2 = false;
    order_r2p1 = false;
    order_r3 = false;
  } else if (order_r_option.compare(ORDER_R_OPTION_R2) == 0) {
    at_least_order_r2 = true;
    order_r2p1 = false;
    order_r3 = false;
  } else if (order_r_option.compare(ORDER_R_OPTION_R2p1) == 0) {
    at_least_order_r2 = true;
    order_r2p1 = true;
    order_r3 = true;
  } else {
    throw std::runtime_error("Invalid setting for order_r_option");
  }
  
}
