#include "discreteconvolution.hpp"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

int main(int /*argc*/, char ** /*argv*/) {

  std::cout << "\nPeriodic convolution: naive vs. fast implementation"
            << std::endl;
  Eigen::VectorXd p(6);
  p << 1, 2, 3, 4, 5, 6;
  Eigen::VectorXd x(6);
  x << 1.1, -1.2, 1.3, -2.4, 2.5, -2.6;
  std::cout << "p\t\t= " << p.transpose() << std::endl;
  std::cout << "x\t\t= " << x.transpose() << std::endl;
  std::cout << "pconv(p,x)\t=" << pconv(p, x).transpose() << std::endl;
  std::cout << "pconv_fast(p,x)\t=" << pconv_fast(p, x).transpose()
            << std::endl;

  return 0;
}
