# include "hermiteapprox.hpp"

int main() {
  auto f = [](double x) { return 1/(1 + x*x); };
  hermiteapprox(f, 0, 1, 50);
  return 0;
}
