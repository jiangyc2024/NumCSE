#ifndef BLOCKDECOMP_HPP
#define BLOCKDECOMP_HPP

#include <iomanip>
#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;

/* \brief Use efficient implementation A*x = bb
 * \param[in] R MatrixXd is nxn and upper triangular
 * \param[in] v VectorXd is nx1
 * \param[in] u VectorXd is nx1
 * \param[in] bb vector is (n+1)x1 and is stacked (b, \beta)^T =: b
 * \param[out] x solution A*bb = x
 */
/* SAM_LISTING_BEGIN_1 */
void solvelse(const MatrixXd& R, const VectorXd& v, const VectorXd& u,
              const VectorXd& bb, VectorXd& x) {
  // Size of R, which is wize of u, v, and size of bb is n+1
  unsigned int n = R.rows();

  // TO DO: (3-9.d)
  // i) Use assert() to check that R, v, u, bb
  // all have the have the appropriate sizes.
  // ii) Use (3-9.b) to solve the LSE.
  // Hint: Use R.triangularView<Upper>() to make
  // use of the triangular structure of R.
  // START

  // END
}
/* SAM_LISTING_END_1 */

/* \brief Use Eigen's LU-solver to solve Ax = y
 * \param[in] R MatrixXd is nxn and upper triangular
 * \param[in] v VectorXd is nx1
 * \param[in] u VectorXd is nx1
 * \param[in] bb vector is (n+1)x1 and is stacked $(b, \beta)^T =: b$
 * \param[out] x solution A*bb = x
 */
/* SAM_LISTING_BEGIN_2 */
bool testSolveLSE(const MatrixXd& R, const VectorXd& v, const VectorXd& u,
                  const VectorXd& b, VectorXd& x) {
  bool areTheSame = false;
  // TO DO: (3-9.e)
  // i) Create the system matrix A and solve the LSE, using
  // an Eigen LU-solver. Store the solution in x.
  // ii) Solve the LSE with solvelse(), and calculate the
  // difference between this solution and x. Return true
  // if and only if the norm of this difference is close
  // enough to zero.
  // START

  // END
  return areTheSame;
}
/* SAM_LISTING_END_2 */

#endif
