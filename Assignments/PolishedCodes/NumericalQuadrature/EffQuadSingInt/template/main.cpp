#include <iostream>

#include "quadsingint.hpp"

///
/// Minimal runner for tasks (8-9.c) and (8-9.e)
///

int main() {
  auto f = [](double t) { return 1. / (2. + std::exp(3 * t)); };

  std::cout << "Integrating 1. / (2. + std::exp(3*t)): " << quadsingint(f, 25)
            << "; should give " << 0.483296828976607 << std::endl;

  tabAndPlotQuadErr();
}
