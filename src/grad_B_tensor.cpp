#include <ctime>
#include <chrono>
#include "qsc.hpp"

#define lp abs_G0_over_B0

using namespace qsc;

/** Compute the grad B tensor
 */
void Qsc::calculate_grad_B_tensor() {
  // In the dimensions of size 3, the order of elements is
  // (normal, binormal, tangent)
  
  qscfloat factor = spsi * B0 / d_l_d_varphi;

  // Evaluate eq (3.12) in Landreman JPP (2021):
  
  // tn
  tempvec = (sG * B0) * curvature;
  grad_B_tensor.set_row(tempvec, 2, 0);

  // nt - same as tn
  grad_B_tensor.set_row(tempvec, 0, 2);

  // bb
  tempvec = factor * (X1c * d_Y1s_d_varphi - iota_N * X1c * Y1c);
  grad_B_tensor.set_row(tempvec, 1, 1);

  // nn
  tempvec = factor * (d_X1c_d_varphi * Y1s + iota_N * X1c * Y1c);
  grad_B_tensor.set_row(tempvec, 0, 0);

  // bn
  tempvec = factor * ((-sG * spsi * d_l_d_varphi) * torsion
		      - iota_N * X1c * X1c);
  grad_B_tensor.set_row(tempvec, 1, 0);

  // nb
  tempvec = factor * (d_Y1c_d_varphi * Y1s - d_Y1s_d_varphi * Y1c
		      + (sG * spsi * d_l_d_varphi) * torsion
		      + iota_N * (Y1s * Y1s + Y1c * Y1c));
  grad_B_tensor.set_row(tempvec, 0, 1);

  // Evaluate eq (3.1) in Landreman JPP (2021):
  for (int j = 0; j < nphi; j++) {
    L_grad_B[j] = B0 * sqrt(2 / (grad_B_tensor(j, 2, 0) * grad_B_tensor(j, 2, 0) +
				 grad_B_tensor(j, 0, 2) * grad_B_tensor(j, 0, 2) +
				 grad_B_tensor(j, 1, 1) * grad_B_tensor(j, 1, 1) +
				 grad_B_tensor(j, 0, 0) * grad_B_tensor(j, 0, 0) +
				 grad_B_tensor(j, 1, 0) * grad_B_tensor(j, 1, 0) +
				 grad_B_tensor(j, 0, 1) * grad_B_tensor(j, 0, 1)));
  }
  L_grad_B_inverse = ((qscfloat)1.0) / L_grad_B;
  grid_min_L_grad_B = L_grad_B.min();
}

void Qsc::calculate_grad_grad_B_tensor() {
  // The elements that follow are computed in the Mathematica notebook "20200407-01 Grad grad B tensor near axis"
  // and then formatted for fortran by the python script process_grad_grad_B_tensor_code

  std::time_t start_time, end_time;
  std::chrono::time_point<std::chrono::steady_clock> start;
  if (verbose > 0) {
    start_time = std::clock();
    start = std::chrono::steady_clock::now();
  }
  if (verbose > 0) std::cout << "Beginning grad_grad_B tensor calculation" << std::endl;
  
  // The order is (normal, binormal, tangent). So element 012 means nbt.

  // Element 000
  work1 = (B0*B0*B0*B0*lp*lp*(8*iota_N*X2c*Y1c*
			      Y1s + 4*iota_N*X2s*
			      (-Y1c*Y1c + Y1s*Y1s) + 
			      2*iota_N*X1c*Y1s*Y20 + 
			      2*iota_N*X1c*Y1s*Y2c - 
			      2*iota_N*X1c*Y1c*Y2s + 
			      5*iota_N*X1c*X1c*Y1c*Y1s*
			      curvature - 
			      2*Y1c*Y20*d_X1c_d_varphi + 
			      2*Y1c*Y2c*d_X1c_d_varphi + 
			      2*Y1s*Y2s*d_X1c_d_varphi + 
			      5*X1c*Y1s*Y1s*curvature*
			      d_X1c_d_varphi + 
			      2*Y1c*Y1c*d_X20_d_varphi + 
			      2*Y1s*Y1s*d_X20_d_varphi - 
			      2*Y1c*Y1c*d_X2c_d_varphi + 
			      2*Y1s*Y1s*d_X2c_d_varphi - 
			      4*Y1c*Y1s*d_X2s_d_varphi))/
    (G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 0, 0, 0);
    
  // Element 001
  work1 = (B0*B0*B0*B0*lp*lp*(Y1c*Y1c*
			      (-6*iota_N*Y2s + 
			       5*iota_N*X1c*Y1s*
			       curvature + 
			       2*(lp*X20*torsion - 
				  lp*X2c*torsion + 
				  d_Y20_d_varphi - 
				  d_Y2c_d_varphi)) + 
			      Y1s*(5*iota_N*X1c*Y1s*Y1s*
				   curvature + 
				   2*(lp*X1c*Y2s*torsion + 
				      Y2s*d_Y1c_d_varphi - 
				      (Y20 + Y2c)*
				      d_Y1s_d_varphi) + 
				   Y1s*(6*iota_N*Y2s + 
					2*lp*X20*torsion + 
					2*lp*X2c*torsion + 
					5*lp*X1c*X1c*curvature*
					torsion + 
					5*X1c*curvature*
					d_Y1c_d_varphi + 
					2*d_Y20_d_varphi + 
					2*d_Y2c_d_varphi)) + 
			      Y1c*(2*(lp*X1c*
				      (-Y20 + Y2c)*torsion - 
				      Y20*d_Y1c_d_varphi + 
				      Y2c*d_Y1c_d_varphi + 
				      Y2s*d_Y1s_d_varphi) + 
				   Y1s*(12*iota_N*Y2c - 
					4*lp*X2s*torsion - 
					5*X1c*curvature*
					d_Y1s_d_varphi - 
					4*d_Y2s_d_varphi))))/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 0, 0, 1);
    
  // Element 002
  work1 = -((B0*B0*B0*lp*lp*(2*Y1c*Y1c*
			     (2*B2c*G0*lp + B0*G2*lp + B0*I2*lp*iota - 
			      2*G0*lp*B20 + 2*B0*G0*iota_N*Z2s + 
			      B0*G0*lp*X20*curvature - 
			      B0*G0*lp*X2c*curvature - 
			      B0*G0*d_Z20_d_varphi + 
			      B0*G0*d_Z2c_d_varphi) + 
			     Y1s*(-2*B0*G0*lp*X1c*Y2s*
				  curvature + 
				  Y1s*(-4*B2c*G0*lp + 2*B0*G2*lp + 
				       2*B0*I2*lp*iota - 4*G0*lp*B20 - 
				       4*B0*G0*iota_N*Z2s + 
				       2*B0*G0*lp*X20*curvature + 
				       2*B0*G0*lp*X2c*curvature + 
				       B0*G0*lp*X1c*X1c*curvature*curvature - 
				       2*B0*G0*d_Z20_d_varphi - 
				       2*B0*G0*d_Z2c_d_varphi)) + 
			     2*G0*Y1c*(B0*lp*X1c*
				       (Y20 - Y2c)*curvature + 
				       2*Y1s*(2*B2s*lp - 2*B0*iota_N*Z2c - 
					      B0*lp*X2s*curvature + 
					      B0*d_Z2s_d_varphi))))/(G0*G0*G0*G0));
  grad_grad_B_tensor.set_row(work1, 0, 0, 2);

  // Element 010
  work1 =-((B0*B0*B0*B0*lp*lp*(3*iota_N*X1c*X1c*X1c*Y1s*
			       curvature + 
			       3*lp*X1c*X1c*Y1s*Y1s*curvature*
			       torsion + 
			       2*(X2s*Y1s*
				  (-2*lp*Y1c*torsion + 
				   d_X1c_d_varphi) + 
				  X20*(lp*Y1c*Y1c*torsion + 
				       lp*Y1s*Y1s*torsion - 
				       Y1c*d_X1c_d_varphi) + 
				  X2c*(-(lp*Y1c*Y1c*
					 torsion) + 
				       lp*Y1s*Y1s*torsion + 
				       Y1c*d_X1c_d_varphi)) - 
			       2*X1c*(3*iota_N*X2s*Y1c - 
				      iota_N*X20*Y1s - 
				      3*iota_N*X2c*Y1s + 
				      lp*Y1c*Y20*torsion - 
				      lp*Y1c*Y2c*torsion - 
				      lp*Y1s*Y2s*torsion - 
				      Y1c*d_X20_d_varphi + 
				      Y1c*d_X2c_d_varphi + 
				      Y1s*d_X2s_d_varphi)))/
	   (G0*G0*G0));
  grad_grad_B_tensor.set_row(work1, 0, 1, 0);

  // Element 011
  work1 =(B0*B0*B0*B0*lp*lp*(-4*iota_N*X1c*Y1s*
			     Y2c + 4*iota_N*X1c*Y1c*
			     Y2s - 3*iota_N*X1c*X1c*Y1c*
			     Y1s*curvature + 
			     2*X20*Y1c*d_Y1c_d_varphi + 
			     2*X20*Y1s*d_Y1s_d_varphi + 
			     3*X1c*X1c*Y1s*curvature*
			     d_Y1s_d_varphi + 
			     2*X2s*(iota_N*Y1c*Y1c - 
				    Y1s*(iota_N*Y1s + 
					 d_Y1c_d_varphi) - 
				    Y1c*d_Y1s_d_varphi) - 
			     2*X2c*(Y1c*
				    (2*iota_N*Y1s + d_Y1c_d_varphi) 
				    - Y1s*d_Y1s_d_varphi) - 
			     2*X1c*Y1c*d_Y20_d_varphi + 
			     2*X1c*Y1c*d_Y2c_d_varphi + 
			     2*X1c*Y1s*d_Y2s_d_varphi))/
    (G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 0, 1, 1);

  // Element 012
  work1 =(2*B0*B0*B0*lp*lp*X1c*
	  (Y1c*(2*B2c*G0*lp + B0*G2*lp + B0*I2*lp*iota - 
		2*G0*lp*B20 + 2*B0*G0*iota_N*Z2s + 
		2*B0*G0*lp*X20*curvature - 
		2*B0*G0*lp*X2c*curvature - 
		B0*G0*d_Z20_d_varphi + 
		B0*G0*d_Z2c_d_varphi) + 
	   G0*Y1s*(2*B2s*lp - 2*B0*iota_N*Z2c - 
		   2*B0*lp*X2s*curvature + 
		   B0*d_Z2s_d_varphi)))/(G0*G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 0, 1, 2);
    
  // Element 020
  work1 =(B0*B0*B0*B0*lp*(-4*lp*lp*X2s*Y1c*Y1s*
			  curvature + 
			  2*lp*lp*X2c*(-Y1c*Y1c + Y1s*Y1s)*
			  curvature + 
			  2*lp*lp*X20*(Y1c*Y1c + Y1s*Y1s)*
			  curvature - 
			  2*lp*lp*X1c*Y1c*Y20*
			  curvature + 
			  2*lp*lp*X1c*Y1c*Y2c*
			  curvature + 
			  2*lp*lp*X1c*Y1s*Y2s*
			  curvature + 
			  3*lp*lp*X1c*X1c*Y1s*Y1s*
			  curvature*curvature + 
			  lp*iota_N*X1c*X1c*X1c*Y1s*
			  torsion - lp*iota_N*X1c*
			  Y1c*Y1c*Y1s*torsion - 
			  lp*iota_N*X1c*Y1s*Y1s*Y1s*
			  torsion - Y1s*Y1s*
			  d_X1c_d_varphi*d_X1c_d_varphi + 
			  iota_N*X1c*X1c*Y1s*
			  d_Y1c_d_varphi - 
			  lp*X1c*Y1s*Y1s*torsion*
			  d_Y1c_d_varphi - 
			  iota_N*X1c*X1c*Y1c*
			  d_Y1s_d_varphi + 
			  lp*X1c*Y1c*Y1s*
			  torsion*d_Y1s_d_varphi + 
			  X1c*Y1s*Y1s*d2_X1c_d_varphi2))/
    (G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 0, 2, 0);

    // Element 021
  work1 =(B0*B0*B0*B0*lp*(-(Y1s*d_X1c_d_varphi*
			    (iota_N*Y1c*Y1c + 
			     Y1s*(iota_N*Y1s + 
				  d_Y1c_d_varphi) - 
			     Y1c*d_Y1s_d_varphi)) + 
			  lp*X1c*X1c*Y1s*
			  (2*iota_N*Y1c*torsion - 
			   torsion*d_Y1s_d_varphi + 
			   Y1s*d_torsion_d_varphi) + 
			  X1c*(Y1c*d_Y1s_d_varphi*
			       (-(iota_N*Y1c) + d_Y1s_d_varphi) 
			       + Y1s*Y1s*(lp*torsion*
					  d_X1c_d_varphi + 
					  iota_N*d_Y1s_d_varphi + 
					  d2_Y1c_d_varphi2) - 
			       Y1s*(d_Y1c_d_varphi*
				    d_Y1s_d_varphi + 
				    Y1c*(-2*iota_N*d_Y1c_d_varphi + 
					 d2_Y1s_d_varphi2)))))/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 0, 2, 1);

  // Element 022
  work1 =(B0*B0*B0*B0*lp*lp*X1c*Y1s*
	  (-(Y1s*curvature*
	     d_X1c_d_varphi) + 
	   X1c*(-(iota_N*Y1c*
		  curvature) + 
		Y1s*d_curvature_d_varphi)))/
    (G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 0, 2, 2);

    // Element 100
  work1 =(-2*B0*B0*B0*B0*lp*lp*X1c*
	  (-2*iota_N*X2s*Y1c + 
	   2*iota_N*X2c*Y1s - 
	   iota_N*X1c*Y2s + 
	   iota_N*X1c*X1c*Y1s*curvature + 
	   lp*X1c*Y1s*Y1s*curvature*
	   torsion - Y20*
	   d_X1c_d_varphi + 
	   Y2c*d_X1c_d_varphi + 
	   Y1c*d_X20_d_varphi - 
	   Y1c*d_X2c_d_varphi - 
	   Y1s*d_X2s_d_varphi))/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 1, 0, 0);

  // Element 101
  work1 =(2*B0*B0*B0*B0*lp*lp*X1c*
	  (lp*X1c*Y20*torsion - 
	   lp*X1c*Y2c*torsion + 
	   Y20*d_Y1c_d_varphi - 
	   Y2c*d_Y1c_d_varphi - 
	   Y2s*d_Y1s_d_varphi + 
	   Y1c*(3*iota_N*Y2s - 
		lp*X20*torsion + 
		lp*X2c*torsion - 
		d_Y20_d_varphi + d_Y2c_d_varphi) 
	   + Y1s*(iota_N*Y20 - 
		  3*iota_N*Y2c - 
		  iota_N*X1c*Y1c*curvature + 
		  lp*X2s*torsion + 
		  X1c*curvature*
		  d_Y1s_d_varphi + d_Y2s_d_varphi))
	  )/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 1, 0, 1);

  // Element 102
  work1 =(2*B0*B0*B0*lp*lp*X1c*
	  (Y1c*(2*B2c*G0*lp + B0*G2*lp + B0*I2*lp*iota - 
		2*G0*lp*B20 + 2*B0*G0*iota_N*Z2s + 
		B0*G0*lp*X20*curvature - 
		B0*G0*lp*X2c*curvature - 
		B0*G0*d_Z20_d_varphi + 
		B0*G0*d_Z2c_d_varphi) + 
	   G0*(B0*lp*X1c*(Y20 - Y2c)*
	       curvature + 
	       Y1s*(2*B2s*lp - 2*B0*iota_N*Z2c - 
		    B0*lp*X2s*curvature + 
		    B0*d_Z2s_d_varphi))))/(G0*G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 1, 0, 2);

  // Element 110
  work1 =(-2*B0*B0*B0*B0*lp*lp*X1c*
	  (lp*X2c*Y1c*torsion + 
	   lp*X2s*Y1s*torsion - 
	   X2c*d_X1c_d_varphi + 
	   X20*(-(lp*Y1c*torsion) + 
		d_X1c_d_varphi) + 
	   X1c*(3*iota_N*X2s + 
		lp*Y20*torsion - 
		lp*Y2c*torsion - 
		d_X20_d_varphi + d_X2c_d_varphi)))/
    (G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 1, 1, 0);

  // Element 111
  work1 =(-2*B0*B0*B0*B0*lp*lp*X1c*
	  (-(iota_N*X2c*Y1s) + 
	   2*iota_N*X1c*Y2s - 
	   X2c*d_Y1c_d_varphi + 
	   X20*(iota_N*Y1s + 
		d_Y1c_d_varphi) + 
	   X2s*(iota_N*Y1c - 
		d_Y1s_d_varphi) - 
	   X1c*d_Y20_d_varphi + 
	   X1c*d_Y2c_d_varphi))/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 1, 1, 1);

  // Element 112
  work1 =(-2*B0*B0*B0*lp*lp*X1c*X1c*
	  (2*B2c*G0*lp + B0*G2*lp + B0*I2*lp*iota - 2*G0*lp*B20 + 
	   2*B0*G0*iota_N*Z2s + 
	   2*B0*G0*lp*X20*curvature - 
	   2*B0*G0*lp*X2c*curvature - 
	   B0*G0*d_Z20_d_varphi + 
	   B0*G0*d_Z2c_d_varphi))/(G0*G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 1, 1, 2);

  // Element 120
  work1 =(B0*B0*B0*B0*lp*X1c*(-2*lp*lp*X20*Y1c*
			      curvature + 
			      2*lp*lp*X2c*Y1c*curvature + 
			      2*lp*lp*X2s*Y1s*curvature + 
			      2*lp*lp*X1c*Y20*curvature - 
			      2*lp*lp*X1c*Y2c*curvature + 
			      2*lp*iota_N*X1c*Y1c*Y1s*
			      torsion - iota_N*X1c*Y1s*
			      d_X1c_d_varphi + 
			      lp*Y1s*Y1s*torsion*
			      d_X1c_d_varphi + 
			      iota_N*X1c*X1c*d_Y1s_d_varphi - 
			      lp*X1c*Y1s*torsion*
			      d_Y1s_d_varphi - 
			      lp*X1c*Y1s*Y1s*
			      d_torsion_d_varphi))/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 1, 2, 0);

  // Element 121
  work1 =(B0*B0*B0*B0*lp*X1c*(-(lp*iota_N*X1c*X1c*
				Y1s*torsion) + 
			      lp*Y1s*torsion*
			      (iota_N*Y1c*Y1c + 
			       Y1s*(iota_N*Y1s + 
				    d_Y1c_d_varphi) - 
			       Y1c*d_Y1s_d_varphi) + 
			      X1c*((iota_N*Y1c - 
				    d_Y1s_d_varphi)*d_Y1s_d_varphi 
				   + Y1s*(-(iota_N*d_Y1c_d_varphi) + 
					  d2_Y1s_d_varphi2))))/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 1, 2, 1);

  // Element 122
  work1 =(B0*B0*B0*B0*lp*lp*X1c*X1c*Y1s*curvature*
	  (iota_N*X1c + 2*lp*Y1s*torsion))/
    (G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 1, 2, 2);

  // Element 200
  work1 =(B0*B0*B0*B0*lp*X1c*Y1s*
	  (lp*iota_N*X1c*X1c*torsion - 
	   lp*iota_N*Y1c*Y1c*torsion - 
	   lp*iota_N*Y1s*Y1s*torsion - 
	   lp*Y1s*torsion*
	   d_Y1c_d_varphi + 
	   X1c*(2*lp*lp*Y1s*curvature*curvature + 
		iota_N*d_Y1c_d_varphi) + 
	   d_X1c_d_varphi*d_Y1s_d_varphi + 
	   Y1c*(iota_N*d_X1c_d_varphi + 
		lp*torsion*d_Y1s_d_varphi) + 
	   Y1s*d2_X1c_d_varphi2))/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 2, 0, 0);

  
  // Element 201
  work1 =(B0*B0*B0*B0*lp*X1c*Y1s*
	  (lp*X1c*(2*iota_N*Y1c*
		   torsion + 
		   Y1s*d_torsion_d_varphi) + 
	   Y1s*(2*lp*torsion*
		d_X1c_d_varphi + 
		2*iota_N*d_Y1s_d_varphi + 
		d2_Y1c_d_varphi2) + 
	   Y1c*(2*iota_N*d_Y1c_d_varphi - 
		d2_Y1s_d_varphi2)))/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 2, 0, 1);

  // Element 202
  work1 =(B0*B0*B0*B0*lp*lp*X1c*X1c*Y1s*
	  (-(iota_N*Y1c*curvature) + 
	   curvature*d_Y1s_d_varphi + 
	   Y1s*d_curvature_d_varphi))/
    (G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 2, 0, 2);

  // Element 210
  work1 =-((B0*B0*B0*B0*lp*X1c*X1c*Y1s*
	    (-2*lp*iota_N*Y1c*torsion + 
	     2*iota_N*d_X1c_d_varphi + 
	     2*lp*torsion*d_Y1s_d_varphi + 
	     lp*Y1s*d_torsion_d_varphi))/
	   (G0*G0*G0));
  grad_grad_B_tensor.set_row(work1, 2, 1, 0);

  // Element 211
  work1 =-((B0*B0*B0*B0*lp*X1c*Y1s*
	    (lp*iota_N*X1c*X1c*torsion - 
	     lp*iota_N*Y1c*Y1c*torsion - 
	     lp*iota_N*Y1s*Y1s*torsion - 
	     lp*Y1s*torsion*
	     d_Y1c_d_varphi - 
	     d_X1c_d_varphi*d_Y1s_d_varphi + 
	     Y1c*(iota_N*d_X1c_d_varphi + 
		  lp*torsion*d_Y1s_d_varphi) + 
	     X1c*(iota_N*d_Y1c_d_varphi - 
		  d2_Y1s_d_varphi2)))/(G0*G0*G0));
  grad_grad_B_tensor.set_row(work1, 2, 1, 1);
  
  // Element 212
  work1 =(B0*B0*B0*B0*lp*lp*X1c*X1c*Y1s*curvature*
	  (iota_N*X1c + 2*lp*Y1s*torsion))/
    (G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 2, 1, 2);

  // Element 220
  work1 =(B0*B0*B0*B0*lp*lp*X1c*X1c*Y1s*
	  (-(iota_N*Y1c*curvature) + 
	   curvature*d_Y1s_d_varphi + 
	   Y1s*d_curvature_d_varphi))/
    (G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 2, 2, 0);

  // Element 221
  work1 =-((B0*B0*B0*B0*lp*lp*X1c*Y1s*curvature*
	    (iota_N*Y1c*Y1c + 
	     Y1s*(iota_N*Y1s + 
		  d_Y1c_d_varphi) - 
	     Y1c*d_Y1s_d_varphi))/(G0*G0*G0));
  grad_grad_B_tensor.set_row(work1, 2, 2, 1);

  // Element 222
  work1 =(-2*B0*B0*B0*B0*lp*lp*lp*X1c*X1c*Y1s*Y1s*
	  curvature*curvature)/(G0*G0*G0);
  grad_grad_B_tensor.set_row(work1, 2, 2, 2);

  // Compute the scale length L_{grad grad B},
  // eq (3.2) in Landreman JPP (2021):
  for (int j = 0; j < nphi; j++) {
    L_grad_grad_B[j] = sqrt(4 * B0 / sqrt(grad_grad_B_tensor(j, 0, 0, 0) * grad_grad_B_tensor(j, 0, 0, 0) +
					  grad_grad_B_tensor(j, 0, 0, 1) * grad_grad_B_tensor(j, 0, 0, 1) +
					  grad_grad_B_tensor(j, 0, 0, 2) * grad_grad_B_tensor(j, 0, 0, 2) +
					  grad_grad_B_tensor(j, 0, 1, 0) * grad_grad_B_tensor(j, 0, 1, 0) +
					  grad_grad_B_tensor(j, 0, 1, 1) * grad_grad_B_tensor(j, 0, 1, 1) +
					  grad_grad_B_tensor(j, 0, 1, 2) * grad_grad_B_tensor(j, 0, 1, 2) +
					  grad_grad_B_tensor(j, 0, 2, 0) * grad_grad_B_tensor(j, 0, 2, 0) +
					  grad_grad_B_tensor(j, 0, 2, 1) * grad_grad_B_tensor(j, 0, 2, 1) +
					  grad_grad_B_tensor(j, 0, 2, 2) * grad_grad_B_tensor(j, 0, 2, 2) +
					  grad_grad_B_tensor(j, 1, 0, 0) * grad_grad_B_tensor(j, 1, 0, 0) +
					  grad_grad_B_tensor(j, 1, 0, 1) * grad_grad_B_tensor(j, 1, 0, 1) +
					  grad_grad_B_tensor(j, 1, 0, 2) * grad_grad_B_tensor(j, 1, 0, 2) +
					  grad_grad_B_tensor(j, 1, 1, 0) * grad_grad_B_tensor(j, 1, 1, 0) +
					  grad_grad_B_tensor(j, 1, 1, 1) * grad_grad_B_tensor(j, 1, 1, 1) +
					  grad_grad_B_tensor(j, 1, 1, 2) * grad_grad_B_tensor(j, 1, 1, 2) +
					  grad_grad_B_tensor(j, 1, 2, 0) * grad_grad_B_tensor(j, 1, 2, 0) +
					  grad_grad_B_tensor(j, 1, 2, 1) * grad_grad_B_tensor(j, 1, 2, 1) +
					  grad_grad_B_tensor(j, 1, 2, 2) * grad_grad_B_tensor(j, 1, 2, 2) +
					  grad_grad_B_tensor(j, 2, 0, 0) * grad_grad_B_tensor(j, 2, 0, 0) +
					  grad_grad_B_tensor(j, 2, 0, 1) * grad_grad_B_tensor(j, 2, 0, 1) +
					  grad_grad_B_tensor(j, 2, 0, 2) * grad_grad_B_tensor(j, 2, 0, 2) +
					  grad_grad_B_tensor(j, 2, 1, 0) * grad_grad_B_tensor(j, 2, 1, 0) +
					  grad_grad_B_tensor(j, 2, 1, 1) * grad_grad_B_tensor(j, 2, 1, 1) +
					  grad_grad_B_tensor(j, 2, 1, 2) * grad_grad_B_tensor(j, 2, 1, 2) +
					  grad_grad_B_tensor(j, 2, 2, 0) * grad_grad_B_tensor(j, 2, 2, 0) +
					  grad_grad_B_tensor(j, 2, 2, 1) * grad_grad_B_tensor(j, 2, 2, 1) +
					  grad_grad_B_tensor(j, 2, 2, 2) * grad_grad_B_tensor(j, 2, 2, 2) ));
  }
  L_grad_grad_B_inverse = ((qscfloat)1.0) / L_grad_grad_B;
  grid_min_L_grad_grad_B = L_grad_grad_B.min();

  if (verbose > 0) {
    end_time = std::clock();
    auto end = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for grad_grad_B from chrono:           "
              << elapsed.count() << std::endl;
    std::cout << "Time for grad_grad_B from ctime (CPU time): "
              << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
  }
}

/** Build the whole grad grad B tensor again using Rogerio's approach,
 * "20200424-01 Rogerio's GradGradB calculation.nb"
 *  This is useful for verifying the two calculations match.
 */
Rank4Tensor Qsc::calculate_grad_grad_B_tensor_alt() {

  Rank4Tensor tensor(nphi, 3, 3, 3);
  
  // The order is (normal, binormal, tangent). So element 012 means nbt.
  
  // Element 000
  work1 =(-2*B0*(-4*sG*spsi*iota_N*X2c*Y1c*	
		 Y1s + iota_N*X1c*X1c*
		 Y1c*(Y1c*
		      (-Y20 + Y2c) + 
		      Y1s*(Y2s - 
			   2*sG*spsi*curvature)) + 
		 X20*Y1c*Y1c*Y1s*
		 d_X1c_d_varphi - 
		 X2c*Y1c*Y1c*Y1s*
		 d_X1c_d_varphi + 
		 X20*Y1s*Y1s*Y1s*
		 d_X1c_d_varphi + 
		 X2c*Y1s*Y1s*Y1s*
		 d_X1c_d_varphi + 
		 sG*spsi*Y1c*Y20*
		 d_X1c_d_varphi - 
		 sG*spsi*Y1c*Y2c*
		 d_X1c_d_varphi - 
		 sG*spsi*Y1s*Y2s*
		 d_X1c_d_varphi - 
		 2*X2s*(sG*spsi*iota_N*Y1s*Y1s + 
			Y1c*Y1c*
			(-(sG*spsi*iota_N) + 
			 iota_N*X1c*Y1s) + 
			Y1c*Y1s*Y1s*
			d_X1c_d_varphi) + 
		 X1c*(iota_N*X2c*Y1c*
		      (-Y1c*Y1c + Y1s*Y1s) + 
		      iota_N*X20*Y1c*
		      (Y1c*Y1c + Y1s*Y1s) - 
		      sG*spsi*iota_N*Y1s*Y20 - 
		      sG*spsi*iota_N*Y1s*Y2c + 
		      sG*spsi*iota_N*Y1c*Y2s - 
		      Y1c*Y1s*Y20*
		      d_X1c_d_varphi + 
		      Y1c*Y1s*Y2c*
		      d_X1c_d_varphi + 
		      Y1s*Y1s*Y2s*
		      d_X1c_d_varphi - 
		      2*sG*spsi*Y1s*Y1s*curvature*
		      d_X1c_d_varphi) - 
		 sG*spsi*Y1c*Y1c*d_X20_d_varphi - 
		 sG*spsi*Y1s*Y1s*d_X20_d_varphi + 
		 sG*spsi*Y1c*Y1c*d_X2c_d_varphi - 
		 sG*spsi*Y1s*Y1s*d_X2c_d_varphi + 
		 2*sG*spsi*Y1c*Y1s*
		 d_X2s_d_varphi))/(lp*spsi);
  tensor.set_row(work1, 0, 0, 0);

  // Element 001
  work1 =(2*B0*(2*iota_N*X2s*Y1c*Y1c*Y1c*
		Y1s + 2*iota_N*X2s*
		Y1c*Y1s*Y1s*Y1s + 
		iota_N*X1c*Y1c*Y1c*Y1c*
		Y20 + iota_N*X1c*Y1c*
		Y1s*Y1s*Y20 - 
		iota_N*X1c*Y1c*Y1c*Y1c*
		Y2c + 6*sG*spsi*iota_N*Y1c*
		Y1s*Y2c - 
		iota_N*X1c*Y1c*Y1s*Y1s*
		Y2c - 3*sG*spsi*iota_N*Y1c*Y1c*
		Y2s - iota_N*X1c*
		Y1c*Y1c*Y1s*Y2s + 
		3*sG*spsi*iota_N*Y1s*Y1s*Y2s - 
		iota_N*X1c*Y1s*Y1s*Y1s*
		Y2s + 2*sG*spsi*iota_N*X1c*
		Y1c*Y1c*Y1s*curvature + 
		2*sG*spsi*iota_N*X1c*Y1s*Y1s*Y1s*
		curvature - 
		2*lp*sG*spsi*X2s*Y1c*
		Y1s*torsion + 
		2*lp*X1c*X2s*Y1c*
		Y1s*Y1s*torsion - 
		lp*sG*spsi*X1c*Y1c*
		Y20*torsion + 
		lp*X1c*X1c*Y1c*Y1s*
		Y20*torsion + 
		lp*sG*spsi*X1c*Y1c*
		Y2c*torsion - 
		lp*X1c*X1c*Y1c*Y1s*
		Y2c*torsion + 
		lp*sG*spsi*X1c*Y1s*
		Y2s*torsion - 
		lp*X1c*X1c*Y1s*Y1s*Y2s*
		torsion + 
		2*lp*sG*spsi*X1c*X1c*Y1s*Y1s*
		curvature*torsion + 
		2*X2s*Y1c*Y1s*Y1s*
		d_Y1c_d_varphi - 
		sG*spsi*Y1c*Y20*
		d_Y1c_d_varphi + 
		X1c*Y1c*Y1s*
		Y20*d_Y1c_d_varphi + 
		sG*spsi*Y1c*Y2c*
		d_Y1c_d_varphi - 
		X1c*Y1c*Y1s*
		Y2c*d_Y1c_d_varphi + 
		sG*spsi*Y1s*Y2s*
		d_Y1c_d_varphi - 
		X1c*Y1s*Y1s*Y2s*
		d_Y1c_d_varphi + 
		2*sG*spsi*X1c*Y1s*Y1s*
		curvature*d_Y1c_d_varphi - 
		2*X2s*Y1c*Y1c*Y1s*
		d_Y1s_d_varphi - 
		X1c*Y1c*Y1c*Y20*
		d_Y1s_d_varphi - 
		sG*spsi*Y1s*Y20*
		d_Y1s_d_varphi + 
		X1c*Y1c*Y1c*Y2c*
		d_Y1s_d_varphi - 
		sG*spsi*Y1s*Y2c*
		d_Y1s_d_varphi + 
		sG*spsi*Y1c*Y2s*
		d_Y1s_d_varphi + 
		X1c*Y1c*Y1s*
		Y2s*d_Y1s_d_varphi - 
		2*sG*spsi*X1c*Y1c*Y1s*
		curvature*d_Y1s_d_varphi + 
		X2c*(Y1c*Y1c - Y1s*Y1s)*
		(iota_N*Y1c*Y1c + iota_N*Y1s*Y1s - 
		 lp*sG*spsi*torsion + 
		 Y1s*(lp*X1c*torsion + 
		      d_Y1c_d_varphi) - 
		 Y1c*d_Y1s_d_varphi) - 
		X20*(Y1c*Y1c + Y1s*Y1s)*
		(iota_N*Y1c*Y1c + iota_N*Y1s*Y1s - 
		 lp*sG*spsi*torsion + 
		 Y1s*(lp*X1c*torsion + 
		      d_Y1c_d_varphi) - 
		 Y1c*d_Y1s_d_varphi) + 
		sG*spsi*Y1c*Y1c*d_Y20_d_varphi + 
		sG*spsi*Y1s*Y1s*d_Y20_d_varphi - 
		sG*spsi*Y1c*Y1c*d_Y2c_d_varphi + 
		sG*spsi*Y1s*Y1s*d_Y2c_d_varphi - 
		2*sG*spsi*Y1c*Y1s*
		d_Y2s_d_varphi))/(lp*spsi);
  tensor.set_row(work1, 0, 0, 1);
  
  // Element 002
  work1 =(-2*(Y1c*Y1c*(G2*spsi + I2*spsi*iota - 
		       2*lp*sG*spsi*B20 + 
		       2*lp*sG*spsi*B2c + 
		       2*B0*sG*spsi*iota_N*Z2s + 
		       B0*lp*sG*spsi*X20*curvature - 
		       B0*lp*sG*spsi*X2c*curvature + 
		       B0*lp*X1c*X20*Y1s*
		       curvature - 
		       B0*lp*X1c*X2c*Y1s*
		       curvature - 
		       B0*sG*spsi*d_Z20_d_varphi + 
		       B0*sG*spsi*d_Z2c_d_varphi) + 
	      Y1s*(B0*lp*X1c*
		   (X20 + X2c)*Y1s*Y1s*
		   curvature - 
		   B0*lp*sG*spsi*X1c*Y2s*
		   curvature + 
		   Y1s*(G2*spsi + I2*spsi*iota - 
			2*lp*sG*spsi*B20 - 
			2*lp*sG*spsi*B2c - 
			2*B0*sG*spsi*iota_N*Z2s + 
			B0*lp*sG*spsi*X20*curvature + 
			B0*lp*sG*spsi*X2c*curvature + 
			B0*lp*X1c*X1c*Y2s*
			curvature + 
			B0*lp*sG*spsi*X1c*X1c*
			curvature*curvature - 
			B0*sG*spsi*d_Z20_d_varphi - 
			B0*sG*spsi*d_Z2c_d_varphi)) + 
	      Y1c*(4*lp*sG*spsi*B2s*
		   Y1s - 
		   B0*(2*lp*X1c*X2s*
		       Y1s*Y1s*curvature + 
		       lp*sG*spsi*X1c*
		       (-Y20 + Y2c)*
		       curvature + 
		       Y1s*
		       (4*sG*spsi*iota_N*Z2c + 
			2*lp*sG*spsi*X2s*curvature + 
			lp*X1c*X1c*Y20*
			curvature - 
			lp*X1c*X1c*Y2c*
			curvature - 
			2*sG*spsi*d_Z2s_d_varphi)))))/
    (lp*spsi);
  tensor.set_row(work1, 0, 0, 2);

  // Element 010
  work1 =(-2*B0*(iota_N*X1c*X1c*X1c*
		 (Y1c*(Y20 - Y2c) + 
		  Y1s*(-Y2s + 
		       sG*spsi*curvature)) - 
		 X1c*X1c*(iota_N*X2c*
			  (-Y1c*Y1c + Y1s*Y1s) + 
			  iota_N*X20*
			  (Y1c*Y1c + Y1s*Y1s) + 
			  Y1s*(-2*iota_N*X2s*
			       Y1c + 
			       lp*(Y1c*
				   (-Y20 + Y2c) + 
				   Y1s*
				   (Y2s - sG*spsi*curvature))*
			       torsion)) + 
		 sG*spsi*(X2s*Y1s*
				  (-2*lp*Y1c*torsion + 
				   d_X1c_d_varphi) + 
				  X20*(lp*Y1c*Y1c*
				       torsion + 
				       lp*Y1s*Y1s*torsion - 
				       Y1c*d_X1c_d_varphi) + 
				  X2c*(-(lp*Y1c*Y1c*
					 torsion) + 
				       lp*Y1s*Y1s*torsion + 
				       Y1c*d_X1c_d_varphi)) + 
		 X1c*(3*sG*spsi*iota_N*X2c*
		      Y1s + 
		      lp*X2c*Y1c*Y1c*Y1s*
		      torsion - 
		      lp*X2c*Y1s*Y1s*Y1s*
		      torsion - 
		      lp*sG*spsi*Y1c*Y20*
		      torsion + 
		      lp*sG*spsi*Y1c*Y2c*
		      torsion + 
		      lp*sG*spsi*Y1s*Y2s*
		      torsion - 
		      X20*Y1s*
		      (-(sG*spsi*iota_N) + 
		       lp*Y1c*Y1c*torsion + 
		       lp*Y1s*Y1s*torsion) + 
		      X2s*Y1c*
		      (-3*sG*spsi*iota_N + 
		       2*lp*Y1s*Y1s*torsion) + 
		      sG*spsi*Y1c*d_X20_d_varphi - 
		      sG*spsi*Y1c*d_X2c_d_varphi - 
		      sG*spsi*Y1s*d_X2s_d_varphi)))/
    (lp*spsi);
  tensor.set_row(work1, 0, 1, 0);
  
  // Element 011
  work1 =(2*B0*(-(X1c*X1c*
		  (Y1c*(Y20 - Y2c) + 
		   Y1s*(-Y2s + 
			sG*spsi*curvature))*
		  (iota_N*Y1c - d_Y1s_d_varphi)) 
		+ X2s*(iota_N*Y1c*Y1c*
		       (sG*spsi - 2*X1c*Y1s) - 
		       sG*spsi*Y1s*
		       (iota_N*Y1s + 
			d_Y1c_d_varphi) + 
		       Y1c*(-(sG*spsi) + 
			    2*X1c*Y1s)*
		       d_Y1s_d_varphi) + 
		sG*spsi*(X20*
				 (Y1c*d_Y1c_d_varphi + 
				  Y1s*d_Y1s_d_varphi) + 
				 X2c*(-(Y1c*
					(2*iota_N*Y1s + 
					 d_Y1c_d_varphi)) + 
				      Y1s*d_Y1s_d_varphi)) + 
		X1c*(-(X2c*
		       (Y1c*Y1c - Y1s*Y1s)*
		       (iota_N*Y1c - 
			d_Y1s_d_varphi)) + 
		     X20*(Y1c*Y1c + Y1s*Y1s)*
		     (iota_N*Y1c - d_Y1s_d_varphi) 
		     + sG*spsi*(Y1c*
					(2*iota_N*Y2s - 
					 d_Y20_d_varphi + 
					 d_Y2c_d_varphi) + 
					Y1s*
					(-2*iota_N*Y2c + 
					 d_Y2s_d_varphi)))))/(lp*spsi);
  tensor.set_row(work1, 0, 1, 1);

  // Element 012
  work1 =(2*X1c*(Y1c*
		 (G2 + I2*iota - 2*lp*sG*B20 + 
		  2*lp*sG*B2c + 
		  2*B0*sG*iota_N*Z2s + 
		  2*B0*lp*sG*X20*curvature - 
		  2*B0*lp*sG*X2c*curvature - 
		  B0*sG*d_Z20_d_varphi + 
		  B0*sG*d_Z2c_d_varphi) + 
		 sG*Y1s*(2*lp*B2s + 
			     B0*(-2*iota_N*Z2c - 
				 2*lp*X2s*curvature + 
				 d_Z2s_d_varphi))))/(lp);
  tensor.set_row(work1, 0, 1, 2);

  // Element 020
  work1 =(B0*(-(lp*sG*spsi*iota_N*Y1c*Y1c*
		torsion) + 
	      lp*iota_N*X1c*X1c*X1c*Y1s*
	      torsion + 
	      X1c*X1c*(lp*lp*Y1s*Y1s*
		       torsion*torsion + 
		       iota_N*Y1s*d_Y1c_d_varphi - 
		       iota_N*Y1c*d_Y1s_d_varphi) + 
	      sG*spsi*Y1c*
	      (iota_N*d_X1c_d_varphi + 
	       2*lp*torsion*d_Y1s_d_varphi) + 
	      X1c*Y1s*
	      (2*lp*lp*sG*spsi*curvature*curvature - 
	       lp*lp*sG*spsi*torsion*torsion - 
	       iota_N*Y1c*d_X1c_d_varphi + 
	       lp*torsion*
	       (Y1s*d_Y1c_d_varphi - 
		Y1c*d_Y1s_d_varphi)) - 
	      Y1s*(Y1s*
		   (lp*sG*spsi*iota_N*torsion + 
		    d_X1c_d_varphi*d_X1c_d_varphi) + 
		   sG*spsi*(2*lp*torsion*
				    d_Y1c_d_varphi - 
				    d2_X1c_d_varphi2))))/(lp*lp*sG);
  tensor.set_row(work1, 0, 2, 0);

  // Element 021
  work1 =(B0*(-(iota_N*Y1c*Y1c*Y1s*
		d_X1c_d_varphi) + 
	      lp*X1c*X1c*Y1s*torsion*
	      (iota_N*Y1c - d_Y1s_d_varphi) + 
	      X1c*(-(iota_N*Y1c*Y1c*
		     d_Y1s_d_varphi) + 
		   Y1c*(lp*sG*spsi*iota_N*
			torsion + 
			iota_N*Y1s*
			d_Y1c_d_varphi + 
			d_Y1s_d_varphi*d_Y1s_d_varphi) - 
		   Y1s*(lp*Y1s*torsion*
			d_X1c_d_varphi + 
			d_Y1c_d_varphi*
			d_Y1s_d_varphi - 
			lp*sG*spsi*d_torsion_d_varphi)) + 
	      Y1s*(-(iota_N*Y1s*Y1s*
		     d_X1c_d_varphi) - 
		   Y1s*d_X1c_d_varphi*
		   d_Y1c_d_varphi + 
		   sG*spsi*(2*lp*torsion*
				    d_X1c_d_varphi + 
				    iota_N*d_Y1s_d_varphi + 
				    d2_Y1c_d_varphi2)) + 
	      Y1c*(sG*spsi*iota_N*
		   d_Y1c_d_varphi + 
		   Y1s*d_X1c_d_varphi*
		   d_Y1s_d_varphi - 
		   sG*spsi*d2_Y1s_d_varphi2)))/
    (lp*lp*sG);
  tensor.set_row(work1, 0, 2, 1);

  // Element 022
  work1 =-((B0*(Y1s*curvature*
		d_X1c_d_varphi + 
		X1c*(iota_N*Y1c*
		     curvature - 
		     Y1s*d_curvature_d_varphi)))/
	   (lp*spsi));
  tensor.set_row(work1, 0, 2, 2);

  // Element 100
  work1 =(-2*B0*X1c*(2*sG*spsi*iota_N*X2c*
		     Y1s + iota_N*X1c*X1c*
		     (Y1c*(Y20 - Y2c) + 
		      sG*spsi*Y1s*curvature) - 
		     X20*Y1c*Y1s*
		     d_X1c_d_varphi + 
		     X2c*Y1c*Y1s*
		     d_X1c_d_varphi - 
		     sG*spsi*Y20*d_X1c_d_varphi + 
		     sG*spsi*Y2c*d_X1c_d_varphi + 
		     X2s*(Y1c*
			  (-2*sG*spsi*iota_N + 
			   iota_N*X1c*Y1s) + 
			  Y1s*Y1s*d_X1c_d_varphi) + 
		     X1c*(-(iota_N*X20*
			    Y1c*Y1c) + 
			  iota_N*X2c*Y1c*Y1c - 
			  sG*spsi*iota_N*Y2s + 
			  lp*sG*spsi*Y1s*Y1s*curvature*
			  torsion + 
			  Y1s*Y20*
			  d_X1c_d_varphi - 
			  Y1s*Y2c*
			  d_X1c_d_varphi) + 
		     sG*spsi*Y1c*d_X20_d_varphi - 
		     sG*spsi*Y1c*d_X2c_d_varphi - 
		     sG*spsi*Y1s*d_X2s_d_varphi))/
    (lp*spsi);
  tensor.set_row(work1, 1, 0, 0);

  // Element 101
  work1 =(-2*B0*X1c*(iota_N*X2s*
		     Y1c*Y1c*Y1s + 
		     iota_N*X2s*Y1s*Y1s*Y1s + 
		     iota_N*X1c*Y1c*Y1c*
		     Y20 - sG*spsi*iota_N*Y1s*
		     Y20 + iota_N*X1c*
		     Y1s*Y1s*Y20 - 
		     iota_N*X1c*Y1c*Y1c*
		     Y2c + 3*sG*spsi*iota_N*Y1s*
		     Y2c - iota_N*X1c*
		     Y1s*Y1s*Y2c - 
		     3*sG*spsi*iota_N*Y1c*Y2s + 
		     sG*spsi*iota_N*X1c*Y1c*
		     Y1s*curvature - 
		     lp*sG*spsi*X2s*Y1s*
		     torsion + 
		     lp*X1c*X2s*Y1s*Y1s*
		     torsion - 
		     lp*sG*spsi*X1c*Y20*
		     torsion + 
		     lp*X1c*X1c*Y1s*Y20*
		     torsion + 
		     lp*sG*spsi*X1c*Y2c*
		     torsion - 
		     lp*X1c*X1c*Y1s*Y2c*
		     torsion + 
		     X2s*Y1s*Y1s*
		     d_Y1c_d_varphi - 
		     sG*spsi*Y20*d_Y1c_d_varphi + 
		     X1c*Y1s*Y20*
		     d_Y1c_d_varphi + 
		     sG*spsi*Y2c*d_Y1c_d_varphi - 
		     X1c*Y1s*Y2c*
		     d_Y1c_d_varphi - 
		     X2s*Y1c*Y1s*
		     d_Y1s_d_varphi - 
		     X1c*Y1c*Y20*
		     d_Y1s_d_varphi + 
		     X1c*Y1c*Y2c*
		     d_Y1s_d_varphi + 
		     sG*spsi*Y2s*d_Y1s_d_varphi - 
		     sG*spsi*X1c*Y1s*
		     curvature*d_Y1s_d_varphi - 
		     X20*Y1c*
		     (iota_N*Y1c*Y1c + iota_N*Y1s*Y1s - 
		      lp*sG*spsi*torsion + 
		      Y1s*(lp*X1c*torsion + 
			   d_Y1c_d_varphi) - 
		      Y1c*d_Y1s_d_varphi) + 
		     X2c*Y1c*
		     (iota_N*Y1c*Y1c + iota_N*Y1s*Y1s - 
		      lp*sG*spsi*torsion + 
		      Y1s*(lp*X1c*torsion + 
			   d_Y1c_d_varphi) - 
		      Y1c*d_Y1s_d_varphi) + 
		     sG*spsi*Y1c*d_Y20_d_varphi - 
		     sG*spsi*Y1c*d_Y2c_d_varphi - 
		     sG*spsi*Y1s*d_Y2s_d_varphi))/
    (lp*spsi);
  tensor.set_row(work1, 1, 0, 1);

  // Element 102
  work1 =(2*X1c*(2*lp*sG*spsi*B2s*
		 Y1s + Y1c*
		 (G2*spsi + I2*spsi*iota - 
		  2*lp*sG*spsi*B20 + 
		  2*lp*sG*spsi*B2c + 
		  2*B0*sG*spsi*iota_N*Z2s + 
		  B0*lp*sG*spsi*X20*curvature - 
		  B0*lp*sG*spsi*X2c*curvature + 
		  B0*lp*X1c*X20*Y1s*
		  curvature - 
		  B0*lp*X1c*X2c*Y1s*
		  curvature - 
		  B0*sG*spsi*d_Z20_d_varphi + 
		  B0*sG*spsi*d_Z2c_d_varphi) - 
		 B0*(lp*X1c*X2s*Y1s*Y1s*
		     curvature + 
		     lp*sG*spsi*X1c*
		     (-Y20 + Y2c)*
		     curvature + 
		     Y1s*(2*sG*spsi*iota_N*Z2c + 
			  lp*sG*spsi*X2s*curvature + 
			  lp*X1c*X1c*Y20*
			  curvature - 
			  lp*X1c*X1c*Y2c*
			  curvature - 
			  sG*spsi*d_Z2s_d_varphi))))/
    (lp*spsi);
  tensor.set_row(work1, 1, 0, 2);

  // Element 110
  work1 =(2*B0*X1c*(iota_N*X1c*X1c*X1c*
		    (Y20 - Y2c) + 
		    X1c*X1c*(-(iota_N*X20*
			       Y1c) + 
			     iota_N*X2c*Y1c + 
			     Y1s*(iota_N*X2s + 
				  lp*(Y20 - Y2c)*
				  torsion)) - 
		    sG*spsi*(lp*X2s*Y1s*
				     torsion + 
				     X2c*(lp*Y1c*
					  torsion - d_X1c_d_varphi) + 
				     X20*(-(lp*Y1c*
					    torsion) + d_X1c_d_varphi)) 
		    + X1c*(-(lp*X20*Y1c*
			     Y1s*torsion) + 
			   lp*X2c*Y1c*Y1s*
			   torsion - 
			   lp*sG*spsi*Y20*torsion + 
			   lp*sG*spsi*Y2c*torsion + 
			   X2s*(-3*sG*spsi*iota_N + 
				lp*Y1s*Y1s*torsion) + 
			   sG*spsi*d_X20_d_varphi - 
			   sG*spsi*d_X2c_d_varphi)))/
    (lp*spsi);
  tensor.set_row(work1, 1, 1, 0);

  // Element 111
  work1 =(-2*B0*X1c*(sG*spsi*
		     (X20 - X2c)*
		     (iota_N*Y1s + d_Y1c_d_varphi) + 
		     X2s*(sG*spsi - 
			  X1c*Y1s)*
		     (iota_N*Y1c - d_Y1s_d_varphi) - 
		     X1c*X1c*(Y20 - Y2c)*
		     (iota_N*Y1c - d_Y1s_d_varphi) + 
		     X1c*(X20*Y1c*
			  (iota_N*Y1c - d_Y1s_d_varphi) 
			  + X2c*Y1c*
			  (-(iota_N*Y1c) + 
			   d_Y1s_d_varphi) + 
			  sG*spsi*(2*iota_N*Y2s - 
					   d_Y20_d_varphi + 
					   d_Y2c_d_varphi))))/(lp*spsi);
  tensor.set_row(work1, 1, 1, 1);

  // Element 112
  work1 =(-2*X1c*X1c*(G2 + I2*iota - 2*lp*sG*B20 + 
		      2*lp*sG*B2c + 2*B0*sG*iota_N*Z2s + 
		      2*B0*lp*sG*X20*curvature - 
		      2*B0*lp*sG*X2c*curvature - 
		      B0*sG*d_Z20_d_varphi + 
		      B0*sG*d_Z2c_d_varphi))/(lp);
  tensor.set_row(work1, 1, 1, 2);

  // Element 120
  work1 =(B0*X1c*(lp*iota_N*Y1c*
		  (sG*spsi + X1c*Y1s)*
		  torsion + 
		  (-(sG*spsi*iota_N) + 
		   lp*Y1s*Y1s*torsion)*
		  d_X1c_d_varphi + 
		  iota_N*X1c*X1c*d_Y1s_d_varphi - 
		  2*lp*sG*spsi*torsion*
		  d_Y1s_d_varphi + 
		  lp*X1c*Y1s*torsion*
		  d_Y1s_d_varphi - 
		  lp*sG*spsi*Y1s*d_torsion_d_varphi))/
    (lp*lp*sG);
  tensor.set_row(work1, 1, 2, 0);

  // Element 121
  work1 =(B0*X1c*(lp*iota_N*Y1c*Y1c*
		  Y1s*torsion + 
		  lp*iota_N*Y1s*Y1s*Y1s*torsion - 
		  lp*lp*sG*spsi*Y1s*torsion*torsion - 
		  sG*spsi*iota_N*d_Y1c_d_varphi + 
		  lp*Y1s*Y1s*torsion*
		  d_Y1c_d_varphi - 
		  lp*Y1c*Y1s*torsion*
		  d_Y1s_d_varphi + 
		  X1c*(-(lp*sG*spsi*iota_N*
			 torsion) + 
		       lp*lp*Y1s*Y1s*torsion*torsion + 
		       (iota_N*Y1c - d_Y1s_d_varphi)*
		       d_Y1s_d_varphi) + 
		  sG*spsi*d2_Y1s_d_varphi2))/(lp*lp*sG);
  tensor.set_row(work1, 1, 2, 1);

  // Element 122
  work1 =(B0*X1c*curvature*
	  (iota_N*X1c + 
	   2*lp*Y1s*torsion))/(lp*spsi);
  tensor.set_row(work1, 1, 2, 2);

  // Element 200
  work1 =(B0*(2*lp*lp*sG*spsi*curvature*curvature + 
	      lp*iota_N*X1c*X1c*torsion - 
	      lp*iota_N*Y1c*Y1c*torsion - 
	      lp*iota_N*Y1s*Y1s*torsion + 
	      iota_N*Y1c*d_X1c_d_varphi + 
	      iota_N*X1c*d_Y1c_d_varphi - 
	      lp*Y1s*torsion*
	      d_Y1c_d_varphi + 
	      lp*Y1c*torsion*
	      d_Y1s_d_varphi + 
	      d_X1c_d_varphi*d_Y1s_d_varphi + 
	      Y1s*d2_X1c_d_varphi2))/
    (lp*lp*spsi);
  tensor.set_row(work1, 2, 0, 0);

  // Element 201
  work1 =(B0*(lp*X1c*(2*iota_N*Y1c*
		      torsion + 
		      Y1s*d_torsion_d_varphi) + 
	      Y1s*(2*lp*torsion*
		   d_X1c_d_varphi + 
		   2*iota_N*d_Y1s_d_varphi + 
		   d2_Y1c_d_varphi2) + 
	      Y1c*(2*iota_N*d_Y1c_d_varphi - 
		   d2_Y1s_d_varphi2)))/(lp*lp*spsi);
  tensor.set_row(work1, 2, 0, 1);

  // Element 202
  work1 =(B0*(-(iota_N*X1c*Y1c*
		curvature) - 
	      Y1s*curvature*
	      d_X1c_d_varphi + 
	      sG*spsi*d_curvature_d_varphi))/(lp*spsi);
  tensor.set_row(work1, 2, 0, 2);

  // Element 210
  work1 =(B0*X1c*(2*lp*iota_N*Y1c*
		  torsion - 
		  2*iota_N*d_X1c_d_varphi - 
		  lp*(2*torsion*d_Y1s_d_varphi + 
		      Y1s*d_torsion_d_varphi)))/
    (lp*lp*spsi);
  tensor.set_row(work1, 2, 1, 0);

  // Element 211
  work1 =-((B0*(iota_N*Y1c*d_X1c_d_varphi + 
		iota_N*X1c*d_Y1c_d_varphi - 
		d_X1c_d_varphi*
		d_Y1s_d_varphi + 
		lp*torsion*
		(iota_N*X1c*X1c - 
		 iota_N*Y1c*Y1c - 
		 Y1s*(iota_N*Y1s + 
		      d_Y1c_d_varphi) + 
		 Y1c*d_Y1s_d_varphi) - 
		X1c*d2_Y1s_d_varphi2))/
	   (lp*lp*spsi));
  tensor.set_row(work1, 2, 1, 1);

  // Element 212
  work1 =(B0*curvature*(iota_N*X1c*X1c + 
			lp*sG*spsi*torsion + 
			lp*X1c*Y1s*torsion))/
    (lp*spsi);
  tensor.set_row(work1, 2, 1, 2);

  // Element 220
  work1 =(B0*(-(iota_N*X1c*Y1c*
		curvature) - 
	      Y1s*curvature*
	      d_X1c_d_varphi + 
	      sG*spsi*d_curvature_d_varphi))/(lp*spsi);
  tensor.set_row(work1, 2, 2, 0);

  // Element 221
  work1 =-((B0*curvature*(iota_N*Y1c*Y1c + 
			  iota_N*Y1s*Y1s - 
			  lp*sG*spsi*torsion + 
			  Y1s*(lp*X1c*torsion + 
			       d_Y1c_d_varphi) - 
			  Y1c*d_Y1s_d_varphi))/
	   (lp*spsi));
  tensor.set_row(work1, 2, 2, 1);

  // Element 222
  work1 =(-2*B0*curvature*curvature)/sG;
  tensor.set_row(work1, 2, 2, 2);

  return tensor;
}
