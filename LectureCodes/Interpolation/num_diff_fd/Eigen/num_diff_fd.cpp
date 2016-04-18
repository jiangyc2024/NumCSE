# include <cmath>
# include "numericaldifferentiation.hpp"

int main() {
  auto exp = [](double x) { return std::exp(x); };
  auto sqrt = [](double x) { return std::sqrt(x); };
  auto dsqrt = [](double x) { return 1./(2*std::sqrt(x)); };
  auto atan = [](double x) { return std::atan(x); };
  auto datan = [](double x) { return 1./(1 + x*x); };

  const double x = 1.1;
  diff(x, exp, exp, "exp");
  diff(x, sqrt, dsqrt, "sqrt");
  diff(x, atan, datan, "atan");

  return 0;
}
