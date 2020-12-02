
#include "gaussConvergence.hpp"

#include <iostream>

int main() {
  // Max. number of nodes
  const unsigned N = 50;

  // "Exact" value of integral
  const double I_ex = 0.870267525725852642;

  // $f(x) = \sinh x$
  std::function<double(double)> f = [](double x) { return std::sinh(x); };

  // PART 1
  double I = gaussConv(f, I_ex, N);

  std::cout << "Approximated integral for " << N
            << "nodes with Gauss-Legendre: " << I;

  // PART 2
  double I_cv = gaussConvCV(f, I_ex, N);

  std::cout << "Approximated integral for " << N
            << "nodes with Gauss-Legendre and variable transformation: "
            << I_cv;

  return 0;
}
