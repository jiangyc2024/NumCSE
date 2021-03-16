#include "quadinf.hpp"

int main() {
  int n = 10;
  double I = quadinf(n, [](double t) { return 1 / (1 + std::pow(t, 2)); });
  // Note: exact integral = pi
  std::cout << "Approximated integral = " << I << std::endl;

  // plot error
  cvgQuadInf();

  return 0;
}
