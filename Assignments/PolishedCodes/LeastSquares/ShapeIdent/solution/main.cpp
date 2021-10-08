#include <Eigen/Dense>
#include <iostream>

#include "shape_ident.hpp"

int main() {
  // Defining shape of stop sign
  Eigen::MatrixXd Xstop(2, 8);
  Xstop << 1, 3, 3, 1, -1, -3, -3, -1, -3, -1, 1, 3, 3, 1, -1, -3;

  // Defining shape of priority sign
  Eigen::MatrixXd Xpriority(2, 8);
  Xpriority << 0, 3, 0, -3, 0, 2.5, 0, -2.5, -3, 0, 3, 0, -2.5, 0, 2.5, 0;

  // First "photo" of a priority sign.
  Eigen::MatrixXd P1(2, 8);
  P1 << 0.23657, 1.35369, -0.13624, -1.33702, 0.0989619, 0.993235, -0.0735973,
      -1.11657, -2.76114, -2.60103, 2.90403, 2.66831, -2.44302, -2.04656,
      2.31922, 2.20296;

  // "Photo" of a sign stop sign.
  Eigen::MatrixXd P2(2, 8);
  P2 << -1.12783, -1.75868, -1.40935, -0.0664574, 1.09654, 1.75873, 1.5195,
      -0.0607661, 1.72169, 0.344036, -0.889686, -1.87847, -1.57535, -0.41511,
      0.834371, 1.88514;

  // "Photo" of a stop sign.
  Eigen::MatrixXd P3(2, 8);
  P3 << -1.23988, -0.731643, 0.00492048, 1.08039, 1.34128, 0.670982, -0.101797,
      -1.02859, 1.53076, 2.02881, 1.36163, -0.340912, -1.47697, -1.99975,
      -1.47947, 0.374859;

  std::cout << "Running tests with:" << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Xstop =" << std::endl;
  std::cout << Xstop << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "XPriority =" << std::endl;
  std::cout << Xpriority << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  Eigen::MatrixXd B = shape_ident_matrix(Xstop);
  std::cout << "Test of shape_ident_matrix(Xstop):" << std::endl;
  std::cout << B << std::endl;
  std::cout << std::endl;

  std::cout << "Testing solve_lsq(Xstop, P1, A):" << std::endl;
  Eigen::MatrixXd A;
  double res = solve_lsq(Xstop, P1, A);
  std::cout << "Norm of the residual: " << res << std::endl;
  std::cout << "A = " << std::endl;
  std::cout << A << std::endl;
  std::cout << std::endl;

  Shape s;
  std::cout << "Testing identify(Xstop, Xpriority, P1, A): " << std::endl;
  s = identify(Xstop, Xpriority, P1, A);
  if (s == Stop) {
    std::cout << "Stop" << std::endl;
  } else {
    std::cout << "Priority" << std::endl;
  }
  std::cout << "Testing identify(Xstop, Xpriority, P2, A):" << std::endl;
  s = identify(Xstop, Xpriority, P2, A);
  if (s == Stop) {
    std::cout << "Stop" << std::endl;
  } else {
    std::cout << "Priority" << std::endl;
  }
  std::cout << "Testing identify(Xstop, Xpriority, P3, A):" << std::endl;
  s = identify(Xstop, Xpriority, P3, A);
  if (s == Stop) {
    std::cout << "Stop" << std::endl;
  } else {
    std::cout << "Priority" << std::endl;
  }

  return 0;
}
