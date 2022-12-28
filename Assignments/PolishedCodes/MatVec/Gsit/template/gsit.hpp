#ifndef GSIT_HPP
#define GSIT_HPP

#include <Eigen/Dense>
#include <algorithm>

/**
 * \brief Use symmetric Gauss-Seidel iterations to solve the system Ax = b
 *
 * \param A system matrix to be decompsed (L + D + U)
 * \param b r.h.s. vector
 * \param x initial guess and last iterate (approximated solution)
 * \param rtol relative tolerance for termination criterion
 */
/* SAM_LISTING_BEGIN_1 */
void GSIt(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
          Eigen::VectorXd& x, double rtol) {
  // TODO: (1-12.c) Implement the Gauss-Seidel iteration.
  // START
  
  // END
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double testGSIt(unsigned int n) {
  double residual_norm = 10.;
  // TODO: (1-12.d) Implement the test from the problem sheet. Return the
  // residual norm.
  // START
  
  // END
  return residual_norm;
}
/* SAM_LISTING_END_2 */

#endif