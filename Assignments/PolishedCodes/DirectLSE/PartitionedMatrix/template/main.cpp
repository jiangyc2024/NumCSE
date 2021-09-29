#include <iostream>

#include "blockdecomp.hpp"

int main() {
  // system dimension
  constexpr unsigned int n = 9;
  // set random seed for reproducibility
  std::srand(9);
  // random test vectors
  Eigen::MatrixXd R =
      Eigen::MatrixXd::Random(n, n).triangularView<Eigen::Upper>();
  Eigen::VectorXd v = Eigen::VectorXd::Random(n);
  Eigen::VectorXd u = Eigen::VectorXd::Random(n);
  Eigen::VectorXd bb = Eigen::VectorXd::Random(n + 1);
  Eigen::VectorXd xe, xo;

  // solve LSE, and test solution
  solvelse(R, v, u, bb, xo);
  const bool works = testSolveLSE(R, v, u, bb, xe);

  std::cout << "LSE solution: " << xo.transpose() << std::endl;

  if (works) {
    std::cout << "Error compared to LU: " << (xe - xo).norm() << std::endl;
  } else {
    std::cout << "testSolveLSE() returns false." << std::endl;
  }
}