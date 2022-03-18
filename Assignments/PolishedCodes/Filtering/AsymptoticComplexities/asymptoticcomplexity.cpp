/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: October 2021
 */

#include <Eigen/Dense>
#include <cassert>
#include <complex>
#include <iostream>
#include <unsupported/Eigen/FFT>

// Some functions for which students should read off the asymptotic complexity.

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd slvtriag(const Eigen::VectorXd &u, const Eigen::VectorXd &v) {
  assert((u.size() == v.size()) && "Size mismatch");
  const Eigen::MatrixXd A{u * v.transpose()};
  return A.triangularView<Eigen::Upper>().solve(u);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd slvsymrom(const Eigen::VectorXd &u) {
  const unsigned int n = u.size();
  return (Eigen::MatrixXd::Identity(n, n) + u * u.transpose()).lu().solve(u);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd getcompbas(const Eigen::MatrixXd &A) {
  const int n = A.rows();
  const int k = A.cols();
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  return qr.householderQ() *
         (Eigen::MatrixXd(n, n - k) << Eigen::MatrixXd::Zero(k, n - k),
          Eigen::MatrixXd::Identity(n - k, n - k))
             .finished();
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXcd fft2_square(const Eigen::MatrixXcd &Y) {
  const unsigned int n = Y.rows();
  assert((n == Y.cols()) && "Matrix Y must be square");
  Eigen::MatrixXcd tmp(n, n);
  Eigen::FFT<double> fft;
  for (int k = 0; k < n; ++k) {
    const Eigen::VectorXcd tv(Y.row(k));
    tmp.row(k) = fft.fwd(tv).transpose();
  }
  for (int k = 0; k < n; ++k) {
    const Eigen::VectorXcd tv(tmp.col(k));
    tmp.col(k) = fft.fwd(tv);
  }
  return tmp;
}
/* SAM_LISTING_END_4 */

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "C++ program for course Numerical Methods for CSE" << std::endl;

  // Call all the above functions
  Eigen::VectorXd u{Eigen::VectorXd::LinSpaced(10, 1.0, 10.0)};
  Eigen::VectorXd v{Eigen::VectorXd::Constant(10, 0.5)};
  std::cout << slvtriag(u, v).transpose() << std::endl;
  std::cout << slvsymrom(u).transpose() << std::endl;
  Eigen::MatrixXd A{Eigen::MatrixXd::Identity(10, 3)};
  std::cout << getcompbas(A) << std::endl;

  return 0;
}
