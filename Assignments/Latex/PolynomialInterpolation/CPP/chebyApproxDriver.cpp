# include "ChebyApprox.hpp"

int main() {
  
  const double maxRho = std::exp(std::asinh(M_PI/3.)), // max rho
               maxRhoAlt = 2.4; // do computations for alternative aswell
  const unsigned N = 1000; // no. of discretization steps

  // Using exact value of maxRho
  std::pair<double, double> res = chebyApprox(maxRho, N);
  std::cout << "using maxRho = " << maxRho << ":\n"
            << "max = " << res.first << ", rho = " << res.second << "\n";

  // Using approximate value of maxRho
  res = chebyApprox(maxRhoAlt, N);
  std::cout << "using maxRho = " << maxRhoAlt << " (suggested alternative):\n"
            << "max = " << res.first << ", rho = " << res.second << "\n";
  return 0;
}
