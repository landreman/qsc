#include <iostream>
#include <ctime>
#include <chrono>
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
  q2.nphi = 31;
  q2.verbose = 0;
  std::time_t start_time, end_time;
  start_time = std::clock();
  auto start = std::chrono::steady_clock::now();
  for (int j = 1; j < 1000; j++) {
    q2.calculate();
  }
  end_time = std::clock();
  auto end = std::chrono::steady_clock::now();

  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Time for 1000 solves from chrono:           "
            << elapsed.count() << std::endl;
  std::cout << "Time for 1000 solves from ctime (CPU time): "
            << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
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
  std::cout << "***********************************************************" << std::endl;
  qsc::Qsc q3;
  q3.input("qsc_in.foo");
  std::cout << "nfp: " << q3.nfp << std::endl;
  std::cout << "eta_bar: " << q3.eta_bar << std::endl;
  q3.calculate();
  std::cout << "Good bye." << std::endl;
  
  return 0;
}
