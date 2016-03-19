# include <iostream>
# include "./diffexcode.hpp"

int main () {
  auto f = [](double x){ return std::exp(x); };
  const double x = -0.5;

  std::cout << "Exact: " << std::exp(x) << "\n"
            << "Approximation: " << diffex(f, x, 0.1, 1e-10, 1e-12) << "\n";

  return 0;
}
