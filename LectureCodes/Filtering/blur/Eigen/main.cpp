#include "blur.hpp"

#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::Index;

int main() {
  const Index n = 11;
  const Index N = 7;
  const MatrixXd P = MatrixXd::Ones(n, n);
  const MatrixXd S = MatrixXd::Ones(N, N);
  const MatrixXd C = blur::blur(P, S);
  std::cout << "C:\n" << C << "\n";
  return 0;
}
