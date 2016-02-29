# include <iostream>
# include <cmath>
# include "simpson.hpp"

int main() {
  auto f = [](double x){ return x*x; };
  auto g = [](double x){ return std::sin(100*x); };
  const double a = 0,
               b = 1;
  const unsigned N = 1000;

  for (unsigned n = 2; n < N; n *= 4){
    std::cout << "=========== " << n << " steps ==========\n";
    std::cout << "x^2 on [0,1]      : " << simpson(f, a, b, n) << "\n";
    std::cout << "sin(100x) on [0,1]: " << simpson(g, a, b, n) << "\n";
  }

  return 0;
}
