#include <iostream>
#include <limits>


#include <Eigen/Dense>
#include "lowrankrep.hpp"

using namespace Eigen;


int main() {
  size_t m = 3;
  size_t n = 2;
  size_t k = 2;
  
  MatrixXd X(m,n);
  X << 5, 0, 2, 1, 7, 4;
  
  std::pair<MatrixXd, MatrixXd> p = factorize_X_AB(X, k);
  
  std::cout << "X = " << std::endl << p.first*p.second.transpose() << std::endl;
  
  MatrixXd A(m,k), B(n,k);
  A << 2, 1, 2, 3, 6, 1;
  B << 4, 4, 5, 0;
  MatrixXd U, S, V;
    
  std::tuple <MatrixXd, MatrixXd, MatrixXd> t = svd_AB(A, B);
    
  std::cout << "U =" << std::endl << std::get<0>(t) << std::endl;
  std::cout << "S =" << std::endl << std::get<1>(t) << std::endl;
  std::cout << "V =" << std::endl << std::get<2>(t) << std::endl;
    
  MatrixXd Ax(m,k), Ay(m,k), Bx(n,k), By(n,k);
  Ax << 1,  0, 9, 2, 6, 3;
  Ay << 8, -2, 3, 4, 5, 8;
  Bx << 2, 1, 2, 3;
  By << 4, 4, 5, 0;
  
  p = rank_k_approx(Ax, Ay, Bx, By);
  
  std::cout << "Z =" << std::endl << p.first*p.second.transpose() << std::endl;
  
}

