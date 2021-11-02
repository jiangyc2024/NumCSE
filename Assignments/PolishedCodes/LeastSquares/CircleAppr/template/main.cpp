///
/// Minimal runner for (9-12)
///

#include <Eigen/Dense>
#include <iostream>

#include "circle_approx.hpp"

int main() {
  Eigen::VectorXd x(8);
  x << 0.7, 3.3, 5.6, 7.5, 6.4, 4.4, 0.3, -1.1;
  Eigen::VectorXd y(8);
  y << 4.0, 4.7, 4.0, 1.3, -1.1, -3.0, -2.5, 1.3;

  compare_convergence(x, y);
  std::cout << "Calculating algebraic fit..." << std::endl;
  const Eigen::Vector3d z_alg = circl_alg_fit(x, y);
  Eigen::Vector3d z_geo_GN = z_alg;
  std::cout << "Calculating geometric fit using Gauss-Newton..." << std::endl;
  circl_geo_fit_GN(x, y, z_geo_GN);
  Eigen::Vector3d z_geo_N = z_alg;
  std::cout << "Calculating geometric fit using Newton..." << std::endl;
  circl_geo_fit_N(x, y, z_geo_N);
  std::cout << "Calculating constrained fit using SVD..." << std::endl;
  const Eigen::Vector3d z_svd = circl_svd_fit(x, y);
  plot(x, y, z_alg, z_geo_GN, z_geo_N, z_svd);
}