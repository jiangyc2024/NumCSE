#include "tridiagqr.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

using namespace Eigen;

int main() {
  std::cout << "\nEnter \"0\" to test all functions.\n"
              << "Enter \"1\" to only test compGivensRotation().\n"
              << "Enter \"2\" to only test TriDiagonalQR::TriDiagonalQR().\n"
              << "Enter \"3\" to only test TriDiagonalQR::applyQT().\n"
              << "Enter \"4\" to only test TriDiagonalQR::solve().\n"
              << "Enter \"5\" to only test invit().\n\n";
  
  int ans=0;
  std::cin >> ans;
  std::cout << std::setprecision(3);
  
  if (ans==1 || ans==0) {
    std::cout << "\e[1mTest for compGivensRotation():\e[0m\n\n";
    MatrixXd a(2,5);
    a << 0, 4, -2,  4, -1.22,
         4, 0,  2, -3, -2.85;
    std::tuple<double,double,double> params;
    Matrix2d G, Grho;
    double rho, gamma, sigma;
    for(int k=0; k<a.cols(); k++) {
      params = compGivensRotation(a.col(k));
      rho = std::get<0>(params);
      gamma = std::get<1>(params);
      sigma = std::get<2>(params);
      G << gamma, sigma, -sigma, gamma;
      Grho = Givens(rho);
      std::cout << "Input vector a = (" << a.col(k).transpose() << ")^T\n";
      std::cout << "Rotation G=\n" << G << std::endl;
      std::cout << "Rotated vector (G^T)*a = (" << a.col(k).transpose()*G << ")^T\n";
      std::cout << "(rho, gamma, sigma) = (" << rho << ", " << gamma << ", " << sigma << ")\n";
      std::cout << "Rotation constructed from rho agrees with G: " << ((G-Grho).norm() < 1E-8) << "\n\n";
    }
    if (ans!=0) return 0;
  }
  int n = 6;
  VectorXd d(n), l(n-1), u(n-1);
  l = VectorXd::LinSpaced(n-1,1,n-1);
  d = VectorXd::LinSpaced(n,n,2*n-1);
  u = VectorXd::LinSpaced(n-1,2*n,3*n-2);
  MatrixXd Adense(n,n);
  Adense.setZero();
  Adense.diagonal() = d;
  Adense.diagonal(-1) = l;
  Adense.diagonal(1) = u;
  std::cout << "Test matrix A=\n" << Adense << "\n\n";
  
  TriDiagonalMatrix A = TriDiagonalMatrix(d,l,u);
  TriDiagonalQR Aqr = TriDiagonalQR(A);
  MatrixXd Q,R;
  auto QR = Aqr.getQRFactors();
	Q = std::get<0>(QR);
	R = std::get<1>(QR);
  
  // Also try a random tridiagonal matrix B.
  std::srand(5);
  l = VectorXd::Random(n-1);
  d = VectorXd::Random(n);
  u = VectorXd::Random(n-1);
  TriDiagonalMatrix B = TriDiagonalMatrix(d,l,u);
  TriDiagonalQR Bqr = TriDiagonalQR(B);
  MatrixXd Bq,Br;
  auto BQR = Bqr.getQRFactors();
	Bq = std::get<0>(BQR);
	Br = std::get<1>(BQR);
  MatrixXd Bdense(n,n);
  Bdense.setZero();
  Bdense.diagonal() = d;
  Bdense.diagonal(-1) = l;
  Bdense.diagonal(1) = u;
  
  if (ans==2 || ans==0) {
    std::cout << "\e[1mTest for TriDiagonalQR::TriDiagonalQR():\e[0m\n\n";
    std::cout << "Computed QR decompostion:\n";
    std::cout << "Q=\n" << Q << std::endl;
    std::cout << "R=\n" << R << std::endl;    
    std::cout << "Q*R=\n" << Q*R << std::endl;
    std::cout << "|QR-A|= " << (Q*R - Adense).norm() << std::endl;
    std::cout << "|Q^TA-R|= " << (Q.transpose()*Adense - R).norm() << "\n\n";
    
    std::cout << "For a random matrix B:\n";
    std::cout << "|QR-B|= " << (Bq*Br - Bdense).norm() << "\n\n";
    if (ans!=0) return 0;
  }
  VectorXd x,y,z;
  if (ans==3 || ans==0) {
    std::cout << "\e[1mTest for TriDiagonalQR::applyQT():\e[0m\n\n";
    
    x = VectorXd::LinSpaced(n,1,n);
    y = Aqr.applyQT(x);
    z = Q.transpose()*x;
    std::cout << "Test vector x = (" << x.transpose() << ")^T\n";
    std::cout << "      (Q^T)*x = (" << y.transpose() << ")^T\n";
    std::cout << "Difference from dense multiplication: " << (z-y).norm() << "\n\n";

    x = VectorXd::Random(n);
    y = Aqr.applyQT(x);
    z = Q.transpose()*x;
    std::cout << "Test vector x = (" << x.transpose() << ")^T\n";
    std::cout << "      (Q^T)*x = (" << y.transpose() << ")^T\n";
    std::cout << "Difference from dense multiplication: " << (z-y).norm() << "\n\n";
    if (ans!=0) return 0;
  }
  if (ans==4 || ans==0) {
    std::cout << "\e[1mTest for TriDiagonalQR::solve():\e[0m\n\n";
    z = VectorXd::LinSpaced(n,1,n);
    
    x = Aqr.solve(z);
    y = Adense.lu().solve(z);
    std::cout << "Right hand side b = (" << z.transpose() << ")^T\n";
    std::cout << "Solution of A*x=b: x = (" << x.transpose() << ")^T\n";
    std::cout << "Difference from dense solver: " << (x-y).norm() << "\n\n";
    
    z = VectorXd::Random(n);
    x = Aqr.solve(z);
    y = Adense.lu().solve(z);
    std::cout << "Right hand side b = (" << z.transpose() << ")^T\n";
    std::cout << "Solution of A*x=b: x = (" << x.transpose() << ")^T\n";
    std::cout << "Difference from dense solver: " << (x-y).norm() << "\n\n";
    
    // We try to solve using a non-invertible matrix B.
    // Set the first line of B to zero.
    std::cout << "Using a non-invertible matrix should throw a runtime error.\n";
    B.d(0) = 0.;
    B.u(0) = 0.;
    Bqr = TriDiagonalQR(B);
    try {
      x = Bqr.solve(z);
    } catch(std::runtime_error& runerr) {
      std::cout << "Non-invertible matrix correctly handled!\n";
      std::cout << "Error message: " << runerr.what() << "\n\n";
    }
    
    if (ans!=0) return 0;
  }
  if (ans==5 || ans==0) {
    std::cout << "\e[1mTest for invit():\e[0m\n\n";
    x = VectorXd::Random(n);
    x *= 10;
    y = x;
    int N, M;
    N = invit<TriDiagonalMatrix>(A, x);
    std::cout << "Using TOL = 1E-6:\n";
    std::cout << "N = " << N << std::endl;
    M = invit<MatrixXd>(Adense, y);
    std::cout << "Difference from dense implementation: " << (N-M) << "\n\n";
    
    x = VectorXd::Random(n);
    x *= 10;
    y = x;
    N = invit<TriDiagonalMatrix>(A, x, 1E-10, 100);
    std::cout << "Using TOL = 1E-10:\n";
    std::cout << "N = " << N << std::endl;
    M = invit<MatrixXd>(Adense, y, 1E-10, 100);
    std::cout << "Difference from dense implementation: " << (N-M) << "\n\n";
    
    if (ans!=0) return 0;
  }
  
  return 0;
}
