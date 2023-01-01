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
  const auto U =
      Eigen::TriangularView<const Eigen::MatrixXd, Eigen::StrictlyUpper>(A);
  const auto L =
      Eigen::TriangularView<const Eigen::MatrixXd, Eigen::StrictlyLower>(A);

  const auto UpD =
      Eigen::TriangularView<const Eigen::MatrixXd, Eigen::Upper>(A);
  const auto LpD =
      Eigen::TriangularView<const Eigen::MatrixXd, Eigen::Lower>(A);

  // A temporary vector to store result of iteration
  Eigen::VectorXd temp(x.size());

  // We'll use pointer magic too
  Eigen::VectorXd* xold = &x;
  Eigen::VectorXd* xnew = &temp;

  double err;

  do {
    // Compute next iteration step
    *xnew = UpD.solve(b) - UpD.solve(L * LpD.solve(b - U * (*xold)));

    // Absolute error
    err = (*xold - *xnew).norm();

    // Swap role of previous/next iteration
    std::swap(xold, xnew);
  } while (err > rtol * (*xnew).norm());

  x = *xnew;
  // END
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double testGSIt(unsigned int n) {
  double residual_norm = 10.;
  // TODO: (1-12.d) Implement the test from the problem sheet. Return the
  // residual norm.
  // START
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
  for (unsigned int i = 0; i < n; ++i) {
    if (i > 0) A(i, i - 1) = 2;
    A(i, i) = 3;
    if (i < n - 1) A(i, i + 1) = 1;
  }
  Eigen::VectorXd b = Eigen::VectorXd::Constant(n, 1);

  Eigen::VectorXd x = b;
  GSIt(A, b, x, 1e-8);
  residual_norm = (A * x - b).norm();
  // END
  return residual_norm;
}
/* SAM_LISTING_END_2 */

#endif