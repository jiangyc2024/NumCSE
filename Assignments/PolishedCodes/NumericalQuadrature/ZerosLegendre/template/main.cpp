#include <iostream>
#include <Eigen/Dense>

#include "legendre.hpp"

///
/// Minimal runner for tasks (8-10.c,d,e,f)
///

int main() {
  constexpr unsigned int n = 8;
  
  // Secant method without regula falsi
  std::cout << "---> Secant method without regula falsi\n";
  Eigen::MatrixXd zeros = gaussPts(n);
  std::cout << "Zeros:\n" << zeros << "\n";
  
  // TODO: (8-10.e) Compute and prin the values of the Legendre polynomials in the zeros obtained.
  // START
  
  // END
  
  // Secant method with regula falsi
  std::cout << "---> Secant method with regula falsi\n";
  
  zeros = gaussPts_regulaFalsi(n);
  std::cout << "Zeros:\n" << zeros << "\n";
  
  for (unsigned int k = 1; k <= n; ++k) {
    Eigen::VectorXd xi = zeros.block(0, k - 1, k, 1);
    Eigen::MatrixXd Lx(k, n + 1), DLx(k, n + 1);
    legvals(xi, Lx, DLx);
    std::cout << "Values of the " << k
    << "-th polynomial in the calculated zeros:\n"
    << Lx.col(k).transpose() << "\n";
  }
}
