#include <Eigen/Dense>
#include <iostream>

#include "matode.hpp"

int main() {
  int n = 3;
  Eigen::MatrixXd A(n, n);
  A << 0, 1, 0, 1, 0, 1, 1, 1, 0;
  Eigen::MatrixXd Y0 = Eigen::MatrixXd::Identity(n, n);

  // test single steps
  std::cout << eeulstep(A, Y0, 0.1) << std::endl
            << std::endl
            << ieulstep(A, Y0, 0.1) << std::endl
            << std::endl
            << impstep(A, Y0, 0.1) << std::endl
            << std::endl;

  std::cout << "Evolution of norm(Y_k'*Y_k - I) for three methods:"
            << std::endl;
  std::tuple<double, double, double> t = checkOrthogonality();

  // Test preservation of orthogonality
  double tol = 1e-13;
  std::cout << "Orthogonality test" << std::endl;
  std::cout << "Explicit Euler: " << (std::get<0>(t) < tol) << std::endl
            << "Implicit Euler: " << (std::get<1>(t) < tol) << std::endl
            << "Implicit Midpoint: " << (std::get<2>(t) < tol) << std::endl;
}
