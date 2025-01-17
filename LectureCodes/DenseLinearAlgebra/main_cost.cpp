// **********************************************************************
// Codes for NumCSE: computational costs of numerical linear algebra
// operations.
// **********************************************************************

#include <Eigen/Dense>
#include <cassert>
#include <iostream>

/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd diagmodsolve1(Eigen::MatrixXd A, const Eigen::VectorXd &b) {
  const Eigen::Index n = A.cols();
  assert((A.rows() == n) && (b.size() == n));
  Eigen::VectorXd x{Eigen::VectorXd::Zero(n)};
  double tmp = A(0, 0);
  for (Eigen::Index i = 0; i < n; ++i) {
    if (i > 0) {
      A(i - 1, i - 1) = tmp;
    }
    tmp = A(i, i);
    A(i, i) *= 2.0;
    x += A.lu().solve(b);
  }
  return x;
}

/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd diagmodsolve2(const Eigen::MatrixXd &A,
                              const Eigen::VectorXd &b) {
  const Eigen::Index n = A.cols();
  assert((A.rows() == n) && (b.size() == n));
  const auto Alu = A.lu(); // \Label[line]{dmslv:1}
  const auto z = Alu.solve(b); 
  const auto W = Alu.solve(Eigen::MatrixXd::Identity(n, n)); // \Label[line]{dmslv:2}
  const Eigen::VectorXd alpha = Eigen::VectorXd::Constant(n, 1.0) +
                                A.diagonal().cwiseProduct(W.diagonal()); 
  if ((alpha.cwiseAbs().array() < 1E-12).any()) {
    throw std::runtime_error("Tiny pivot!");
  }
  return n * z - W * z.cwiseProduct(A.diagonal().cwiseQuotient(alpha));
}

/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double sumtrv1(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
  const Eigen::Index n = A.cols();
  assert((A.rows() == n) && (b.size() == n));

  return b.transpose() *
         A.triangularView<Eigen::Upper>().solve(
             Eigen::MatrixXd::Identity(n, n)) *
         b;
}

/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
double sumtrv2(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
  const Eigen::Index n = A.cols();
  assert((A.rows() == n) && (b.size() == n));

  return b.transpose() * A.triangularView<Eigen::Upper>().solve(b);
}

/* SAM_LISTING_END_3 */

int main() {
  const Eigen::MatrixXd A{
      (Eigen::MatrixXd(3, 3) << 3., 2., 1., 5., 6., 4., 7., 8., 9.).finished()};
  const Eigen::VectorXd b = (Eigen::VectorXd(3) << 1., 2., 3.).finished();
  std::cout << "DMS1 = " << diagmodsolve1(A, b).transpose() << std::endl;
  std::cout << "DMS2 = " << diagmodsolve2(A, b).transpose() << std::endl;
  std::cout << "TI1 = " << sumtrv1(A, b) << std::endl;
  std::cout << "TI2 = " << sumtrv2(A, b) << std::endl;
  return 0;
}
