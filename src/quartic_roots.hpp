#include "vector_matrix.hpp"

/** Find the roots of a quartic equation, using the companion matrix.
 *
 *  This subroutine follows the same algorithm as Matlab's 'roots' routine.
 *  coefficients should have 5 elements. They are ordered the same way as in matlab, from the coefficient of x^4 to the coefficient of x^0.
 *  real_parts and imag_parts should have 4 elements.
 */
void quartic_roots(qsc::qscfloat* coefficients, qsc::qscfloat* real_parts, qsc::qscfloat* imag_parts);
