# include <iostream>
# include "./splineapprox.hpp"

int main() {
  auto f = [](double x) { return std::atan(x); };
  auto df = [](double x) { return 1/(1 + x*x); };

  auto g = [](double x) { return 1/(1 + std::exp(-2*x)); };
  auto dg = [](double x) { return 2*exp(-2*x)/std::pow(1 + exp(-2*x), 2); }; 

  std::cout << "Spline approximation for f = atan: \n";
  splineapprox(f, df, -5, 5, 100, "f");

  std::cout << "Spline approximation for f = 1/(1 + exp(-2x)): \n";
  splineapprox(g, dg, -1, 1, 100, "g");

  return 0;
}
