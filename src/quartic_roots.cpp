#include "quartic_roots.hpp"

using namespace qsc;

// Representation of the LAPACK routine we need:
extern "C" {
  void dgeev_(char* JOBVL, char* JOBVR, int* N, double* A, int* LDA,
	      double* WR, double* WI, double* VL, int* LDVL, double* VR,
	      int* LDVR, double* WORK, int* IWORK, int* INFO);
  
  void sgeev_(char* JOBVL, char* JOBVR, int* N, float* A, int* LDA,
	      float* WR, float* WI, float* VL, int* LDVL, float* VR,
	      int* LDVR, float* WORK, int* IWORK, int* INFO);
}

// Choose either the single or double precision version of BLAS/LAPACK routines:
#ifdef SINGLE

#define geev_ sgeev_

#else

#define geev_ dgeev_

#endif

/** Find the roots of a quartic equation, using the companion matrix.
 *
 *  coefficients should have 5 elements.
 *  real_parts and imag_parts should have 4 elements.
 */
void quartic_roots(qscfloat* coefficients, qscfloat* real_parts, qscfloat* imag_parts) {

  qscfloat matrix[16];
  const int LWORK = 100; // Size of work array for LAPACK
  int IWORK = LWORK;
  qscfloat WORK[LWORK], VL[1], VR[1];
  int INFO, four = 4, one = 1;
  char JOBVR = 'N';

  matrix[0]  = -coefficients[1] / coefficients[0];
  matrix[1]  = 1.0;
  matrix[2]  = 0.0;
  matrix[3]  = 0.0;

  matrix[4]  = -coefficients[2] / coefficients[0];
  matrix[5]  = 0.0;
  matrix[6]  = 1.0;
  matrix[7]  = 0.0;
  
  matrix[8]  = -coefficients[3] / coefficients[0];
  matrix[9]  = 0.0;
  matrix[10] = 0.0;
  matrix[11] = 1.0;

  matrix[12] = -coefficients[4] / coefficients[0];
  matrix[13] = 0.0;
  matrix[14] = 0.0;
  matrix[15] = 0.0;

  geev_(&JOBVR, &JOBVR, &four, matrix, &four, real_parts, imag_parts,
	VL, &one, VR, &one, WORK, &IWORK, &INFO);
  
  if (INFO != 0) {
    std::cerr << "Error in DGEEV: info=" << INFO << std::endl;
    std::cerr << "coefficients:" << coefficients[0] << " " << coefficients[1] << " " << coefficients[2] << " " << coefficients[3] << " " << coefficients[4] << std::endl;
    std::cerr << "Top row of companion matrix:" << matrix[0] << " " << matrix[4] << " " << matrix[8] << " " << matrix[12] << std::endl;
    throw std::runtime_error("LAPACK error in *geev");
  }
}

  /*
*> \brief <b> DGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGEEV + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeev.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeev.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeev.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
*                         LDVR, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   WI( * ), WORK( * ), WR( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEEV computes for an N-by-N real nonsymmetric matrix A, the
*> eigenvalues and, optionally, the left and/or right eigenvectors.
*>
*> The right eigenvector v(j) of A satisfies
*>                  A * v(j) = lambda(j) * v(j)
*> where lambda(j) is its eigenvalue.
*> The left eigenvector u(j) of A satisfies
*>               u(j)**H * A = lambda(j) * u(j)**H
*> where u(j)**H denotes the conjugate-transpose of u(j).
*>
*> The computed eigenvectors are normalized to have Euclidean norm
*> equal to 1 and largest component real.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBVL
*> \verbatim
*>          JOBVL is CHARACTER*1
*>          = 'N': left eigenvectors of A are not computed;
*>          = 'V': left eigenvectors of A are computed.
*> \endverbatim
*>
*> \param[in] JOBVR
*> \verbatim
*>          JOBVR is CHARACTER*1
*>          = 'N': right eigenvectors of A are not computed;
*>          = 'V': right eigenvectors of A are computed.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the N-by-N matrix A.
*>          On exit, A has been overwritten.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] WR
*> \verbatim
*>          WR is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] WI
*> \verbatim
*>          WI is DOUBLE PRECISION array, dimension (N)
*>          WR and WI contain the real and imaginary parts,
*>          respectively, of the computed eigenvalues.  Complex
*>          conjugate pairs of eigenvalues appear consecutively
*>          with the eigenvalue having the positive imaginary part
*>          first.
*> \endverbatim
*>
*> \param[out] VL
*> \verbatim
*>          VL is DOUBLE PRECISION array, dimension (LDVL,N)
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored one
*>          after another in the columns of VL, in the same order
*>          as their eigenvalues.
*>          If JOBVL = 'N', VL is not referenced.
*>          If the j-th eigenvalue is real, then u(j) = VL(:,j),
*>          the j-th column of VL.
*>          If the j-th and (j+1)-st eigenvalues form a complex
*>          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
*>          u(j+1) = VL(:,j) - i*VL(:,j+1).
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the array VL.  LDVL >= 1; if
*>          JOBVL = 'V', LDVL >= N.
*> \endverbatim
*>
*> \param[out] VR
*> \verbatim
*>          VR is DOUBLE PRECISION array, dimension (LDVR,N)
*>          If JOBVR = 'V', the right eigenvectors v(j) are stored one
*>          after another in the columns of VR, in the same order
*>          as their eigenvalues.
*>          If JOBVR = 'N', VR is not referenced.
*>          If the j-th eigenvalue is real, then v(j) = VR(:,j),
*>          the j-th column of VR.
*>          If the j-th and (j+1)-st eigenvalues form a complex
*>          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
*>          v(j+1) = VR(:,j) - i*VR(:,j+1).
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR.  LDVR >= 1; if
*>          JOBVR = 'V', LDVR >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,3*N), and
*>          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
*>          performance, LWORK must generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  if INFO = i, the QR algorithm failed to compute all the
*>                eigenvalues, and no eigenvectors have been computed;
*>                elements i+1:N of WR and WI contain eigenvalues which
*>                have converged.
		    */
