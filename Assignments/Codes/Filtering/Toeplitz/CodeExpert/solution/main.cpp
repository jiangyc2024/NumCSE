#include <Eigen/Dense>
#include <iostream>

#include "timer.h"
#include "toeplitz.hpp"

using namespace Eigen;

int main() {
  int n;

  // Initialization
  n = 3;
  VectorXd c(n), r(n), x(n);
  c << 1, 2, 3;
  r << 1, 4, 5;
  x << 6, 7, 8;

  // Compute with both functions toepmatmult and toepmult
  std::cout << "Check that toepmatmult and toepmult are correct" << std::endl;
  VectorXd y_1 = toepmatmult(c, r, x);
  VectorXd y_2 = toepmult(c, r, x);
  std::cout << "Error = " << (y_1 - y_2).norm() << std::endl;

  // Initialization
  n = 4;
  VectorXd h(n), y(n);
  h << 1, 2, 3, 4;
  y << 5, 6, 7, 8;

  // Compute with both functions ttmatsolve and ttrecsolve
  std::cout << "Check that ttmatsolve and ttrecsolve are correct" << std::endl;
  VectorXd x_1 = ttmatsolve(h, y);
  VectorXd x_2 = ttrecsolve(h, y, 2);
  std::cout << "Error = " << (x_1 - x_2).norm() << std::endl;

  // compute and print the runtime comparsion for toepmatmult vs toepmult 
  // and for ttmatsolve vs ttrecsole
  std::cout << "Runtime comparison of "
            << "toepmatmult vs toepmult and ttmatsolve vs ttrecsolve"
            << std::endl;


  // call runtime computation
  runtime_toeplitz();

}
