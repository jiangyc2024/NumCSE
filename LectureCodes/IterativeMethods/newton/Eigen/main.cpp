#include "newton.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

/* SAM_LISTING_BEGIN_0 */
void newton2Ddriver(void) {
  // Function F defined through lambda function
  auto F = [](const Eigen::Vector2d &x) {
    Eigen::Vector2d z;
    const double x1 = x(0), x2 = x(1);
    z << x1 * x1 - 2 * x1 - x2 + 1, x1 * x1 + x2 * x2 - 1;
    return (z);
  };
  // Lambda function for the computation of the Newton correction
  auto DFinv = [](const Eigen::Vector2d &x, const Eigen::Vector2d &f) {
    Eigen::Matrix2d J;
    const double x1 = x(0), x2 = x(1);
    // Jacobian of F
    J << 2 * x1 - 2, -1, 2 * x1, 2 * x2;
    // Solve 2x2 linear system of equations
    Eigen::Vector2d s = J.lu().solve(f);
    return (s);
  };
  // initial guess
  Eigen::Vector2d x0(2., 3.);

  // Invoke Newton's method
  Eigen::Vector2d x = newton(F, DFinv, x0, 1E-6, 1E-8);
  std::cout << "||F(x)|| = " << F(x).norm() << std::endl;
}
/* SAM_LISTING_END_0 */

int main() { newton2Ddriver(); }
