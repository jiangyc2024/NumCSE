# include <cmath>
# include "numericaldifferentiation.hpp"

// Numerical differentiation for \Blue{$f(x) = e^x$}, \Blue{$g(x) = \sqrt{x}$} and \Blue{$h(x) = \arctan{x}$}
int main() {
  auto exp = [](double x) { return std::exp(x); };
  auto sqrt = [](double x) { return std::sqrt(x); };
  auto dsqrt = [](double x) { return 1./(2*std::sqrt(x)); };
  auto atan = [](double x) { return std::atan(x); };
  auto datan = [](double x) { return 1./(1 + x*x); };

  const double x = 1.1;
  diff(x, exp, exp, "exp"); // d(exp) = exp
  diff(x, sqrt, dsqrt, "sqrt");
  diff(x, atan, datan, "atan");

  return 0;
}
