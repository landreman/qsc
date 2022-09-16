#include <cmath>
#include "qsc.hpp"

using namespace qsc;

void qsc::Qsc::mercier() {

  qscfloat abs_G0 = std::abs(G0);
  
  // See Overleaf note "Mercier criterion near the magnetic axis- detailed notes".
  // See also "20200604-02 Checking sign in Mercier DGeod near axis.docx"

  // integrand = d_l_d_phi0 * (Y1c * Y1c + X1c * (X1c + Y1s)) / (Y1c * Y1c + (X1c + Y1s) * (X1c + Y1s))
  work1 = d_l_d_phi0 * (eta_bar*eta_bar*eta_bar*eta_bar + curvature*curvature*curvature*curvature*sigma*sigma + eta_bar*eta_bar*curvature*curvature) 
    / (eta_bar*eta_bar*eta_bar*eta_bar + curvature*curvature*curvature*curvature*(1+sigma*sigma) + 2*eta_bar*eta_bar*curvature*curvature);

  qscfloat integral = work1.sum() * d_phi0 * nfp * 2 * pi / axis_length;

  // DGeod_times_r2 = -(2 * sign_G * sign_psi * mu0 * mu0 * p2 * p2 * G0 * G0 * G0 * G0 * eta_bar * eta_bar &
  DGeod_times_r2 = -(2 * mu0 * mu0 * p2 * p2 * G0 * G0 * G0 * G0 * eta_bar * eta_bar 
		     / (pi * pi * pi * B0 * B0 * B0 * B0 * B0 * B0 * B0 * B0 * B0 * B0 * iota_N * iota_N))
    * integral;

  d2_volume_d_psi2 = 4*pi*pi*abs_G0/(B0*B0*B0)*(3*eta_bar*eta_bar - 4*B20_mean/B0 + 2*(G2+iota*I2)/G0);
  /*
  std::cout << "etabar term: " << 4*pi*pi*abs_G0/(B0*B0*B0)*(3*eta_bar*eta_bar)
	    << "  B20 term: " << 4*pi*pi*abs_G0/(B0*B0*B0)*( - 4*B20_mean/B0)
	    << "  p2 term: " << 4*pi*pi*abs_G0/(B0*B0*B0)*( 2*(G2+iota*I2)/G0) << std::endl;
  */
  DWell_times_r2 = (mu0 * p2 * abs_G0 / (8 * pi * pi * pi * pi * B0 * B0 * B0)) * (d2_volume_d_psi2 - 8 * pi * pi * mu0 * p2 * abs_G0 / (B0 * B0 * B0 * B0 * B0));

  DMerc_times_r2 = DWell_times_r2 + DGeod_times_r2;

}
