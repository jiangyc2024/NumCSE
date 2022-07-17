#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "matmatCOO.hpp"
#include "matplotlibcpp.h"
#include "timer.h"

namespace plt = matplotlibcpp;

using Trip = Eigen::Triplet<double>;
using TripVec = std::vector<Trip>;

int main() {
  constexpr unsigned int n = 6;
  Eigen::MatrixXd A(n, n), B(n, n);
  A << 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
  B << 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  // COO format
  TripVec A_COO = Mat2COO(A);
  TripVec B_COO = Mat2COO(B);

  // Compute with standard matrix multiplication and both new multipliers
  std::cout << "--> Check that the multipliers are correct" << std::endl;
  TripVec C_eigen, C_naive, C_effic;

  C_eigen = Mat2COO(A * B);
  C_naive = COOprod_naive(A_COO, B_COO);
  C_effic = COOprod_effic(A_COO, B_COO);

  Eigen::MatrixXd C_mat_eigen;
  Eigen::MatrixXd C_mat_naive;
  Eigen::MatrixXd C_mat_effic;
  C_mat_eigen = COO2Mat(C_eigen);
  C_mat_naive = COO2Mat(C_naive);
  C_mat_effic = COO2Mat(C_effic);

  std::cout << "Error eigen vs naive = " << (C_mat_eigen - C_mat_naive).norm()
            << std::endl;
  std::cout << "Error naive vs effic = " << (C_mat_naive - C_mat_effic).norm()
            << std::endl;

  constexpr unsigned int repeats = 3;

  // Compute runtimes of different multipliers for products between sparse
  // matrices
  std::cout << "--> Runtime comparison of naive vs efficient multiplier"
            << std::endl;
  std::cout << "--> Product between sparse matrices" << std::endl;

  std::cout << std::setw(20) << "n" << std::setw(20) << "time naive [s]"
            << std::setw(20) << "time effic [s]" << std::endl;

  std::vector<double> time_naive, time_effic, linear, quadratic, cubic;
  std::vector<double> sizes;

  // Loop over matrix size
  for (unsigned int k = 4; k <= 10; ++k) {
    Timer tm_naive, tm_effic;
    const unsigned int n = std::pow(2, k);
    sizes.push_back(n);

    // Repeat test multiple times
    for (unsigned int r = 0; r < repeats; ++r) {
      // Initialization of random sparse matrices
      A = randMat(n, n, 0.1);
      B = randMat(n, n, 0.1);

      // COO format
      TripVec A_COO = Mat2COO(A);
      TripVec B_COO = Mat2COO(B);

      // Compute runtime with naive solver
      tm_naive.start();
      C_naive = COOprod_naive(A_COO, B_COO);
      tm_naive.stop();
      // Compute runtime with efficient solver
      tm_effic.start();
      C_effic = COOprod_effic(A_COO, B_COO);
      tm_effic.stop();
    }

    linear.push_back(n);
    quadratic.push_back(n * n);
    cubic.push_back(n * n * n);
    time_naive.push_back(tm_naive.min());
    time_effic.push_back(tm_effic.min());

    // Print runtimes
    std::cout << std::setw(20) << n << std::scientific << std::setprecision(3)
              << std::setw(20) << tm_naive.min() << std::setw(20)
              << tm_effic.min() << std::endl;
  }

  plt::figure();
  plt::loglog(sizes, linear, "k:", {{"label", "O(n)"}});
  plt::loglog(sizes, quadratic, "k-", {{"label", "O(n^2)"}});
  plt::loglog(sizes, cubic, "k--", {{"label", "O(n^3)"}});
  plt::loglog(sizes, time_naive, "r+", {{"label", "naive"}});
  plt::loglog(sizes, time_effic, "b+", {{"label", "efficient"}});
  plt::xlabel("Matrix size (n)");
  plt::ylabel("Time [s]");
  plt::legend();
  plt::title("Comparison of timings - sparse matrices");
  plt::savefig("./cx_out/matmatCOO_comparison_sparse.png");

  sizes.clear();
  linear.clear();
  quadratic.clear();
  cubic.clear();
  time_naive.clear();
  time_effic.clear();

  // Compute runtimes of different multipliers for products between sparse OR
  // dense matrices
  std::cout << "--> Runtime comparison of naive vs efficient multiplier"
            << std::endl;
  std::cout << "--> Product between sparse OR dense matrices" << std::endl;

  std::cout << std::setw(20) << "n" << std::setw(20) << "time naive [s]"
            << std::setw(20) << "time effic [s]" << std::endl;

  // Loop over matrix size
  for (unsigned int k = 4; k <= 8; ++k) {
    Timer tm_naive, tm_effic;
    const unsigned int n = std::pow(2, k);
    sizes.push_back(n);

    // Repeat test multiple times
    for (unsigned int r = 0; r < repeats; ++r) {
      // Initialization of random sparse matrices
      A = Eigen::MatrixXd::Random(n, n);
      B = Eigen::MatrixXd::Random(n, n);

      // COO format
      TripVec A_COO = Mat2COO(A);
      TripVec B_COO = Mat2COO(B);

      // Compute runtime with naive solver
      tm_naive.start();
      C_naive = COOprod_naive(A_COO, B_COO);
      tm_naive.stop();
      // Compute runtime with efficient solver
      tm_effic.start();
      C_effic = COOprod_effic(A_COO, B_COO);
      tm_effic.stop();
    }

    linear.push_back(n);
    quadratic.push_back(n * n);
    cubic.push_back(n * n * n);
    time_naive.push_back(tm_naive.min());
    time_effic.push_back(tm_effic.min());

    // Print runtimes
    std::cout << std::setw(20) << n << std::scientific << std::setprecision(3)
              << std::setw(20) << tm_naive.min() << std::setw(20)
              << tm_effic.min() << std::endl;
  }

  plt::figure();
  plt::loglog(sizes, linear, "k:", {{"label", "O(n)"}});
  plt::loglog(sizes, quadratic, "k-", {{"label", "O(n^2)"}});
  plt::loglog(sizes, cubic, "k--", {{"label", "O(n^3)"}});
  plt::loglog(sizes, time_naive, "r+", {{"label", "naive"}});
  plt::loglog(sizes, time_effic, "b+", {{"label", "efficient"}});
  plt::xlabel("Matrix size (n)");
  plt::ylabel("Time [s]");
  plt::legend();
  plt::title("Comparison of timings - sparse/dense matrices");
  plt::savefig("./cx_out/matmatCOO_comparison_sparsedense.png");
}
