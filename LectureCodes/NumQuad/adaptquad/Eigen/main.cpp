# include "./adaptquad.hpp"
# include <Eigen/Dense>
# include <cmath>
# include <iostream>

using Eigen::VectorXd;

int main () {
  auto f =[](double x) { return std::exp(-x*x); };
  VectorXd M(4);
  M << -100, 0.1, 0.5, 100;
  std::cout << "Sqrt(Pi) - Int_{-100}^{100} exp(-x*x) dx = ";
  std::cout << adaptquad::adaptquad(f, M, 1e-10, 1e-12) - std::sqrt(M_PI) << "\n";
  return 0;
}
