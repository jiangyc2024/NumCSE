#include <Eigen/Dense>
#include <iostream>

#include "rosenbrock.hpp"

int main() {
  auto f = [](const Eigen::Vector3d y) {
    Eigen::Vector3d fy;
    fy << y(0) * y(1), y(1) * y(2), y(2) - y(0);
    return fy;
  };
  auto df = [](const Eigen::Vector3d y) {
    Eigen::Matrix3d J;
    J << y(1), y(0), 0, 0, y(2), y(1), -1, 0, 1;
    return J;
  };
  Eigen::Vector3d y0;
  y0 << 1, 2, 3;
  std::cout << "Test of solveRosenbrock():"
            << solveRosenbrock(f, df, y0, 10, 2.).at(10).transpose()
            << std::endl;

  double cvgRate = cvgRosenbrock();
  std::cout << "Convergence rate: " << std::round(cvgRate) << std::endl;
}
