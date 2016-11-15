//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

// Evaluate the Legendre polynomials and its derivatives in vector $x$ 
// using the 3-term recursion formulae.
// The outputs are the matrices $Lx$ and $DLx$.
void legvals(const VectorXd& x, MatrixXd& Lx, MatrixXd& DLx) {
  const int n = Lx.cols()-1;
  const int N = x.size();
    // TODO: evaluate Legendre polynomials
}

// Evaluate $P_n(x)$ for a scalar $x$ and integer $n$.
double Pnx(const double x, const int n) {
  VectorXd Px(n+1);
    // TODO: evaluate $P_n(x)$
  return Px(n);
}

// Find the Gauss points using the secant method without regula falsi.
MatrixXd gaussPts(const int n, const double rtol=1e-10,
							   const double atol=1e-12) {
  MatrixXd zeros(n,n);
    // TODO: secant method without regula falsi
  return zeros;
}

// Find the Gauss points using the secant method with regula falsi.
// The standard secant method may be obtained
// by commenting out line 106.
MatrixXd gaussPts_regulaFalsi(const int n, const double rtol=1e-10,
										   const double atol=1e-12) {
  MatrixXd zeros(n,n);
    // TODO: secant method with regula falsi
  return zeros;
}

int main() {
  const int n = 8;
  MatrixXd zeros;
  
  // Secant method without regula falsi
  std::cout << "---> Secant method without regula falsi\n";
  
  zeros = gaussPts(n);
  std::cout << "Zeros:\n" << zeros << "\n";
    
  for (int k = 1; k < n+1; ++k) {
    VectorXd xi = zeros.block(0, k-1, k, 1);
    MatrixXd Lx(k,n+1), DLx(k,n+1);
    legvals(xi, Lx, DLx);
    std::cout << "Values of the " << k
			  << "-th polynomial in the calculated zeros:\n"
              << Lx.col(k).transpose() << "\n";
  }
  
  // Secant method with regula falsi
  std::cout << "---> Secant method with regula falsi\n";
  
  zeros = gaussPts_regulaFalsi(n);
  std::cout << "Zeros:\n" << zeros << "\n";
    
  for (int k = 1; k < n+1; ++k) {
    VectorXd xi = zeros.block(0, k-1, k, 1);
    MatrixXd Lx(k,n+1), DLx(k,n+1);
    legvals(xi, Lx, DLx);
    std::cout << "Values of the " << k
			  << "-th polynomial in the calculated zeros:\n"
              << Lx.col(k).transpose() << "\n";
  }
}
