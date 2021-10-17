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
  u << 1, 2, 3, 4;
  v << 10, 20, 30;

  // Compute with both functions
  std::cout << "Check that all functions are correct" << std::endl;

  Eigen::VectorXd uv_1 = polyMult_naive(u, v);
  std::cout << "Naive multiplicator: " << std::endl << uv_1 << std::endl;

  Eigen::VectorXd uv_2 = polyMult_fast(u, v);

  // polyMult_fast() returns an uninitialized vector if not implemented i.e. a
  // vector of size 0. Trying to print a vector of size 0 results in an error.
  //
  if (uv_2.size() == 0) {
    std::cout << "polyMult_fast() not implemented yet.";
  } else {
    std::cout << "Efficient multiplicator: " << std::endl << uv_2 << std::endl;
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
    for (int p = 2; p <= 15; ++p) {
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

  return 0;
}
