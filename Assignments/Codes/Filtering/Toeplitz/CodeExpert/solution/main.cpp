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

  // Initialization
  int repeats = 3;
  VectorXd out;

  // Compute runtimes of different multiplicators
  std::cout << "Runtime comparison of "
            << "toepmatmult vs toepmult and ttmatsolve vs ttrecsolve"
            << std::endl;

  // Loop over vector size
  for (int l = 3; l <= 11; ++l) {
    // Timers
    Timer tm_matmult, tm_mult, tm_ttmat, tm_ttrec;
    int n = pow(2, l);

    // Repeat test many times
    for (int repeat = 0; repeat < repeats; ++repeat) {
      c = VectorXd::Random(n);
      r = VectorXd::Random(n);
      r(0) = c(0);
      x = VectorXd::Random(n);
      h = VectorXd::LinSpaced(n, 1, n).cwiseInverse();
      y = VectorXd::Random(n);

      // Compute runtime of toepmatmult
      tm_matmult.start();
      out = toepmatmult(c, r, x);
      tm_matmult.stop();
      // Compute runtime of toepmult
      tm_mult.start();
      out = toepmult(c, r, x);
      tm_mult.stop();
      // Compute runtime of ttmatsolve
      tm_ttmat.start();
      out = ttmatsolve(h, y);
      tm_ttmat.stop();
      // Compute runtime of ttrecsolve
      tm_ttrec.start();
      out = ttrecsolve(h, y, l);
      tm_ttrec.stop();
    }

    // Print progress
    std::cout << n << " completed" << std::endl;
  }
}
