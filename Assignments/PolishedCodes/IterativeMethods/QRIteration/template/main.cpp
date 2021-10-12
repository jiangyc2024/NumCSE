#include "qriteration.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <vector>

int main(int /*argc*/, char ** /*argv*/) {

  {
    std::cout << "Givens rotation for [3,4]" << std::endl;
    Eigen::Vector2d a(3.0, 4.0);
    Eigen::Vector2d g = givens(a);
    Eigen::Matrix2d G(2, 2);
    G << g[0], -g[1], g[1], g[0];
    a.applyOnTheLeft(G);
    std::cout << "G*a = " << a.transpose() << std::endl << std::endl;
  }

  unsigned n = 5;
  Eigen::VectorXd dT = Eigen::VectorXd::LinSpaced(n, 1, n);
  Eigen::VectorXd uT = 0.1 * Eigen::VectorXd::LinSpaced(n - 1, 1, n - 1);
  Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n - 1; ++i) {
    T(i, i) = dT[i];
    T(i, i + 1) = T(i + 1, i) = uT[i];
  }
  T(n - 1, n - 1) = dT[n - 1];
  std::cout << "Symmetric tridiagonal matrix T:= \n"
            << T << std::endl;

  Eigen::VectorXd dTp, uTp;
  std::tie(dTp, uTp) = qrStep(dT, uT);
  Eigen::MatrixXd Tp = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n - 1; ++i) {
    Tp(i, i) = dTp[i];
    Tp(i, i + 1) = Tp(i + 1, i) = uTp[i];
  }
  Tp(n - 1, n - 1) = dTp[n - 1];
  std::cout << "Obtained symmetric tridiagonal T' from qrStep(d(T), u(T)): \n"
            << Tp << std::endl;

  Eigen::HouseholderQR<Eigen::MatrixXd> qr(T);
  const Eigen::MatrixXd Rfull =
      qr.matrixQR().template triangularView<Eigen::Upper>();

  const ::Eigen::MatrixXd Qfull =
      qr.householderQ() * Eigen::MatrixXd::Identity(n, n);
  const ::Eigen::MatrixXd Tpfull = Rfull * Qfull;

  // NOTE: The matrices obtained from qrStep() and QR decomposition may agree
  //       only in absolute values.
  std::cout << "Obtained symmetric tridiagonal T' from QR decomposition: \n"
            << Tpfull << std::endl;

  return 0;
}
