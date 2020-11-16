#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "adaptivepolyintp.hpp"

int main() {
  auto f = [](double t) { return std::sin(std::exp(2 * t)); };

  // Test interval
  const double a = 0, b = 1;

  // Get interpolation nodes and print runtimes
  const unsigned N = 100;   // no. of sampling points
  const double tol = 1e-4;  // tolerance

  Eigen::VectorXd tf;  // nodes for interpolation of f
  tf = adaptivepolyintp(f, a, b, tol, N);

  std::cout << tf << std::endl << std::endl;

  std::vector<double> ef;  // errors for f
  tf = adaptivepolyintp(f, a, b, tol, N, &ef);

  std::cout << "n \t\t \\eps_n" << std::endl
            << "-------------------------" << std::endl;
  for (unsigned int i = 0; i < ef.size(); ++i) {
    std::cout << i + 1 << "\t\t" << ef[i] << std::endl;
  }

  plotInterpolationError();
  return 0;
}
