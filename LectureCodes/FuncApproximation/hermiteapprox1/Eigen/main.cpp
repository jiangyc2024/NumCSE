# include "hermiteapprox1.hpp"
# include <cmath>

int main() {
  auto f = [](double x) { return 1/(1 + x*x); };
  auto df = [](double x) { return -1/std::pow(1 + x*x, 2)*2*x; };
  hermiteapprox(f, df, 0, 1, 50);
  return 0;
}
