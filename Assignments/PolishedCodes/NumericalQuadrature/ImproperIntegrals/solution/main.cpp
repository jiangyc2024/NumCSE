#include "quadinf.hpp"

///
/// Minimal runner for both tasks (8-3.e, 8-3.f)
///

int main() {
  constexpr int n = 10;
  double I = quadinf(n, [](double t) { return 1 / (1 + std::pow(t, 2)); });
  // Note: exact integral = pi
  std::cout << "Approximated integral = " << I << std::endl;

  // plot error
  cvgQuadInf();

  return 0;
}
