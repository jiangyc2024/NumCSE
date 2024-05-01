#include <cmath>
#include <functional>
#include <iostream>

#include "odeintssctrl.h"

int main() {
  // A scalar initial-value problem
  using State_t = double;
  using DiscEvolOp = std::function<State_t(double, State_t)>;
  // Definition of right-hand-side for ODE
  const auto f = [](double y) { return y * y; };
  // Initial value
  const double y0 = 0.5;
  // Exact solution
  auto y = [](double t) { return 1.0 / (2.0 - t); };
  // The norm is just the modulus in the scalar case
  const auto norm = [](double x) { return std::fabs(x); };

  // Low-order method: explicit euler (order 1)
  const DiscEvolOp psilow = [&](double h, double y) { return y + h * f(y); };

  // "High-order method": explicit trapezoidal rule (order 2)
  const DiscEvolOp psihigh = [&](double h, double y) {
    const double k1 = f(y);
    const double k2 = f(y + h * k1);
    return y + (h / 2.) * (k1 + k2);
  };
  // Call simple adaptive integrator
  const std::vector<std::pair<double, double>> states =
      odeintssctrl(psilow, 1, psihigh, y0, {.T=1.9, .h0=0.2, .reltol=1e-3, .abstol=1e-4, .hmin=1e-4}, norm);
  // Output solution and error
  std::cout << "Adaptive integration, " << states.size() - 1 << " timesteps"
            << std::endl;
  for (const auto &ty : states) {
    std::cout << "t = " << ty.first << ": y = " << ty.second
              << ", error = " << norm(ty.second - y(ty.first)) << std::endl;
  }
  return 0;
}
