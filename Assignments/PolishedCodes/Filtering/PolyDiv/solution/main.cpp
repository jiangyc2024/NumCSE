#include <iomanip>
#include <iostream>

#include "polydiv.hpp"
#include "timer.h"

int main() {
  // Initialization
  int m = 4;
  int n = 3;
  Eigen::VectorXd u(m);
  Eigen::VectorXd v(n);
  u << 1, 2, 3, 4;  // u = 1 + 2x + 3x^2 + 4x^3
  v << 10, 20, 30;  // v = 10 + 20x + 30x^2

  // Compute with both functions
  std::cout << "Check that all functions are correct" << std::endl;

  Eigen::VectorXd uv_1 = polyMult_naive(u, v);
  std::cout << "Naive multiplicator: " << std::endl << uv_1 << std::endl;

  Eigen::VectorXd uv_2 = polyMult_fast(u, v);
  std::cout << "Efficient multiplicator: " << std::endl << uv_2 << std::endl;

  if (uv_1.size() != uv_2.size()) {
    std::cout << "Vectors returned from polyMult_naive() and polyMult_fast()"
              << "have different sizes: " << uv_1.size() << " vs. "
              << uv_2.size() << std::endl
              << "Fix it." << std::endl;
    return 0;
  }

  std::cout << "Error = " << (uv_1 - uv_2).norm() << std::endl;

  Eigen::VectorXd v_new = polyDiv(uv_2, u);
  std::cout << "Error of efficient division = " << (v - v_new).norm()
            << std::endl;

  // Initialization
  int repeats = 3;
  Eigen::VectorXd uv;

  // Compute runtimes of different multiplicators
  std::cout << "Runtime comparison of naive v efficient multiplicator"
            << std::endl;

  // Header
  std::cout << std::setw(20) << "n" << std::setw(20) << "time naive [s]"
            << std::setw(20) << "time fast [s]" << std::endl;

  // Loop over vector size
  for (int p = 2; p <= 10; ++p) {
    // Timers
    Timer tm_naive, tm_effic;
    int n = pow(2, p);

    // Repeat test many times
    for (int r = 0; r < repeats; ++r) {
      u = Eigen::VectorXd::Random(n);
      v = Eigen::VectorXd::Random(n);

      // Compute runtime of naive multiplicator
      tm_naive.start();
      uv = polyMult_naive(u, v);
      tm_naive.stop();
      // Compute runtime of efficient multiplicator
      tm_effic.start();
      uv = polyMult_fast(u, v);
      tm_effic.stop();
    }

    // Print runtimes
    std::cout << std::setw(20) << n << std::scientific << std::setprecision(3)
              << std::setw(20) << tm_naive.min() << std::setw(20)
              << tm_effic.min() << std::endl;
  }
}
