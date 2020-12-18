#include <iostream>
#include "qsc.hpp"

int main(int argc, char* argv[]) {
  qsc::Matrix m = qsc::differentiation_matrix(4, 0, 2 * qsc::pi);
  std::cout << m << std::endl;

  m = 4.3;
  std::cout << m(1, 1) << std::endl;
  std::cout << m << std::endl;

  /*
  qsc::Qsc q;
  q.nfp = 3;
  q.nphi = 63;
  
  q.R0c.resize(2, 0.0);
  q.R0s.resize(2, 0.0);
  q.Z0c.resize(2, 0.0);
  q.Z0s.resize(2, 0.0);
  q.R0c[0] = 1.0;
  q.R0c[1] = 0.045;
  q.Z0s[1] = -0.045;
  
  q.allocate();
  q.init_axis();
  */
  qsc::Qsc q("r1 section 5.1");
  //std::cout << "d_d_phi:" << std::endl << q.d_d_phi;
  //std::cout << "d_d_varphi:" << std::endl << q.d_d_varphi;
  std::cout << "***********************************************************" << std::endl;
  qsc::Qsc q2("1");
  //std::cout << "d_d_phi:" << std::endl << q2.d_d_phi;
  //std::cout << "d_d_varphi:" << std::endl << q2.d_d_varphi;
  /*
  qsc::Vector foo;
  foo.resize(3, 42.0);
  */
  
  //qsc::Matrix foo(3, 3);
  //  qsc::Vector foo(6);
  //  std::vector<double> foo(5);
  //std::cout << foo << std::endl;
  /*
  std::cout << "R0:" << q.R0 << std::endl;
  std::cout << "Z0:" << q.Z0 << std::endl;
  */
  std::cout << "Good bye." << std::endl;
  
  return 0;
}
