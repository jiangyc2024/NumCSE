#include <Eigen/Dense>
#include <iostream>

#include "lowrankrep.hpp"

int main() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  constexpr std::size_t k = 2;

  Eigen::MatrixXd X(3, 2);
  X << 5, 0, 2, 1, 7, 4;

  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> p = factorize_X_AB(X, k);

  std::cout << "X = " << std::endl
            << p.first * p.second.transpose() << std::endl;

  Eigen::MatrixXd A(m, k), B(n, k);
  A << 2, 1, 2, 3, 6, 1;
  B << 4, 4, 5, 0, 7, 4, 3, 2;
  Eigen::MatrixXd U, S, V;

  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> t =
      svd_AB(A, B);

  std::cout << "U =" << std::endl << std::get<0>(t) << std::endl;
  std::cout << "S =" << std::endl << std::get<1>(t) << std::endl;
  std::cout << "V =" << std::endl << std::get<2>(t) << std::endl;

  Eigen::MatrixXd Ax(6, k), Ay(6, k), Bx(5, k), By(5, k);
  Ax << 0, 9, 2, 6, 3, 5, -2, 3, 4, 8, 9, 0;
  Ay << 8, -2, 3, 4, 5, 8, 7, 5, 6, 3, 2, 1;
  Bx << 2, 1, 2, -3, 6, 7, 8, 4, 5, 1;
  By << 4, 4, -5, 0, 3, 2, 5, 9, 0, 5;

  p = rank_k_approx(Ax, Ay, Bx, By);

  std::cout << "Z =" << std::endl
            << p.first * p.second.transpose() << std::endl;
}