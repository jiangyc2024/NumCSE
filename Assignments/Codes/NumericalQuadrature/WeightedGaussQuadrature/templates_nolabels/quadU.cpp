//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>

// Function implementing weighted Gauss quadrature
template<typename Function>
double quadU(const Function& f, const unsigned n) {
  double q = 0, w, xi;
  for (unsigned j = 0; j < n; j++) {
    w = M_PI/(n+1) * std::pow(std::sin((j+1)*M_PI/(n+1)), 2);
    xi = std::cos( (j+1.)/(n+1)*M_PI );
    q += w*f(xi);
  }
  return q;
}

// Test the implementation. The parameter q of the exponential decay is approximated by 0.1-0.2. 
// After the 18th iteration only numerical error is present.
int main(){
  auto f = [](double t) { return 1. / (2+std::exp(3*t)); };
  const double exact = 0.483296828976607;
  std::vector<double> e(25);
  std::cout << std::setw(5) << "n" << std::setw(20) << "Error" << std::setw(20) << "Approximated q" << "\n";
  for (unsigned int n = 0; n < 25; n++) {
    e[n] = std::abs(quadU(f,n+1)-exact);
    if (n > 1)
     std::cout << std::setw(5) << n << std::setw(20) << e[n] << std::setw(20) << e[n]/e[n-1] << "\n";
  }

  return 0;
}
