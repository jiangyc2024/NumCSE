
#include "adaptedlinreg.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
using namespace Eigen;

int main() {
  int n = 25;
  VectorXd t, x, y, noise;
  t = VectorXd::LinSpaced(n,0.0,1.0);
  std::srand(41);
  noise = VectorXd::Random(n);
  
  std::cout << "Linear problem, y = alpha+beta*t\n";
  y = 12*t + VectorXd::Constant(n,-154) + noise;
  x = linReg(t, y);
  std::cout << std::setprecision(6);
  std::cout << "(alpha, beta) = " << x.transpose() << std::endl;
  
  std::cout << "\nNonlinear problem, y = alpha*exp(beta*t)\n";
  y = 17*exp( -3*t.array() ).matrix() + 0.5*noise;
  x = expFit(t, y);
  std::cout << "(alpha, beta) = " << x.transpose() << std::endl;
}
