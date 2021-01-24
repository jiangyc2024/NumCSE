#include "gaussConvergence.hpp"

#include <iostream>
#include <iomanip>

int main() {
  // Max. number of nodes
  constexpr unsigned N = 50;

  // "Exact" value of integral
  constexpr double I_ex = 0.870267525725852642;

  // $f(x) = \sinh x$
  std::function<double(double)> f = [](double x) { return std::sinh(x); };
  
  std::cout << "Exact value of integral: " << std::setprecision(10) << I_ex << std::endl;

  // PART 1
  const double I = gaussConv(f, I_ex, N);

  std::cout << "Approximated integral for " << N
  << " nodes with Gauss-Legendre: " << std::setprecision(10) << I << std::endl;

  // PART 2
  const double I_cv = gaussConvCV(f, I_ex, N);

  std::cout << "Approximated integral for " << N
            << " nodes with Gauss-Legendre and variable transformation: "
            << std::setprecision(10) << I_cv << std::endl;

  return 0;
}
