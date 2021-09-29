#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "matplotlibcpp.h"
#include "structuredLSE.hpp"
#include "timer.h"

namespace plt = matplotlibcpp;

int main() {
  unsigned int n = 9;
  // Compute with both solvers
  std::cout << "--> Check that the solvers are correct" << std::endl;
  Eigen::VectorXd a = Eigen::VectorXd::Random(n);
  Eigen::VectorXd b = Eigen::VectorXd::Random(n);
  Eigen::VectorXd x1, x2, x;
  std::vector<double> naive(6), fast(6), on3(6), on(6), n_(6);

  solveA(a, b, x1);
  solveA_fast(a, b, x2);

  std::cout << "Error = " << (x1 - x2).norm() << std::endl;

  // Compute runtimes of different solvers
  std::cout << "--> Runtime comparison of naive vs fast solver" << std::endl;
  // Number of repetitions
  constexpr unsigned int repeats = 3;

  // Header
  std::cout << std::setw(20) << "n" << std::setw(20) << "time naive [s]"
            << std::setw(20) << "time fast [s]" << std::endl;

  unsigned int i = 0;
  // Loop over matrix size
  for (n = 16; n <= 512; n <<= 1, ++i) {
    // Timers
    Timer tm_naive, tm_fast;

    // Repeat test many times
    for (unsigned int r = 0; r < repeats; ++r) {
      a = Eigen::VectorXd::Random(n);
      b = Eigen::VectorXd::Random(n);

      // Compute runtime with naive solver
      tm_naive.start();
      solveA(a, b, x);
      tm_naive.stop();
      // Compute runtime with efficient solver
      tm_fast.start();
      solveA_fast(a, b, x);
      tm_fast.stop();
    }

    naive[i] = tm_naive.min();
    fast[i] = tm_fast.min();
    on3[i] = 1e-7 * std::pow(n, 3);
    on[i] = 1e-7 * n;
    n_[i] = n;

    // Print runtimes
    std::cout << std::setw(20) << n << std::scientific << std::setprecision(3)
              << std::setw(20) << naive[i] << std::setw(20) << fast[i]
              << std::endl;
  }

  // naive plot
  plt::figure();
  plt::loglog(n_, naive, "r+", {{"label", "original"}});
  plt::loglog(n_, on3, "--", {{"label", "O(n^3)"}});
  plt::xlabel("Vector size (n)");
  plt::ylabel("Time [s]");
  plt::legend();
  plt::title("Timings of naive solver");
  plt::savefig("cx_out/naive.png");

  // comparison plot
  plt::figure();
  plt::loglog(n_, naive, "r+", {{"label", "original"}});
  plt::loglog(n_, fast, "b+", {{"label", "efficient"}});
  plt::loglog(n_, on3, "--", {{"label", "O(n^3)"}});
  plt::loglog(n_, on, "-", {{"label", "O(n)"}});
  plt::xlabel("Vector size (n)");
  plt::ylabel("Time [s]");
  plt::legend();
  plt::title("Comparison of timings");
  plt::savefig("cx_out/comparison.png");
}