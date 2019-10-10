
#include "tridiagleastsquares.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace Eigen;

int main() {
  std::cout << "Enter problem dimension (n>=4): ";
  int n;
  std::cin >> n;
  // Define tridiagonal matrix.
  double alpha = -4.0;
  double beta = 1.0;
  SparseMatrix<double> T(n,n);
  std::vector<Triplet<double>> triplets;
  triplets.reserve(3*n-2);
  for(int i=0; i<(n-1); i++) {
    triplets.push_back(Triplet<double>(i,i,alpha));
    triplets.push_back(Triplet<double>(i+1,i,beta));
    triplets.push_back(Triplet<double>(i,i+1,beta));
  }
  triplets.push_back(Triplet<double>(n-1,n-1,alpha));
  T.setFromTriplets(triplets.begin(), triplets.end());
  T.makeCompressed();
  
  // Synthetic data
  VectorXd z = VectorXd::LinSpaced(n,-1.0,1.0);
  std::srand(39);
  VectorXd c = T*z + 0.5*VectorXd::Random(n);
  VectorXd x = lsqEst(z,c);
  
  std::cout << "Least squares parameters: \n"
            << "(alpha, beta) = " << std::setprecision(6)
            << x.transpose() << std::endl;
}