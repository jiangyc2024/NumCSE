#include <Eigen/Dense>
#include <iostream>

#include "ArrowMatrix.hpp"

int main(void) {
  // Test vectors
  Eigen::VectorXd a(5);
  a << 1., 2., 3., 4., 5.;
  Eigen::VectorXd d(5);
  d << 1., 3., 4., 5., 6.;
  Eigen::VectorXd x(5);
  x << -5., 4., 6., -8., 5.;
  Eigen::VectorXd yi;

  // Run both functions
  arrow_matrix_2_times_x(a, d, x, yi);
  Eigen::VectorXd ye(yi.size());
  efficient_arrow_matrix_2_times_x(a, d, x, ye);

  // Compute error
  const double err = (yi - ye).norm();

  // Output error
  std::cout << "--> Correctness test." << std::endl;
  std::cout << "Error: " << err << std::endl;

  // Print out runtime
  std::cout << "--> Runtime test." << std::endl;
  runtime_arrow_matrix();
  std::cout << "Plot was created. See 'Files'." << std::endl;
}
