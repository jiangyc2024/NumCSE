///
/// Minimal runner for (8-8)
///

#include <iostream>

#include "julia.hpp"

int main() {
  // Tests will only pass if z = (-0.888,1.333).
  Eigen::Vector2d z(-0.888, 1.333);
  std::cout << "Test: F(z) = " << juliaF(z).transpose() << "\n\n";
  std::cout << "Test: DF(z) = \n" << juliaDF(z) << "\n\n";

  // No tests.
  julia();
}
