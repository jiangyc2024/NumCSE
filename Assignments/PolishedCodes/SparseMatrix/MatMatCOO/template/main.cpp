#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "matmatCOO.hpp"
#include "matplotlibcpp.h"
#include "timer.h"

namespace plt = matplotlibcpp;

using Trip = Eigen::Triplet<double>;
using TripVec = std::vector<Trip>;

/**
 * @brief Build matrix $A$ as Eigen::MatrixXd from COO format
 *
 * @param A An $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @return Eigen::MatrixXd The $m \times n$ matrix from vector of triplets $A$
 */
Eigen::MatrixXd COO2Mat(const TripVec& A) {
  unsigned int m = 0, n = 0;
  for (auto const& a : A) {
    if (a.row() > m) {
      m = a.row();
    }
    if (a.col() > n) {
      n = a.col();
    }
  }
  ++m, ++n;  // first index is 0
  Eigen::MatrixXd A_mat = Eigen::MatrixXd::Zero(m, n);

  for (auto const& a : A) {
    A_mat(a.row(), a.col()) += a.value();
  }

  return A_mat;
}

/**
 * @brief Build random binary matrix $A$ as Eigen::MatrixXd
 *
 * @param m Number of desired rows for $A$
 * @param n Number of desired cols for $A$
 * @param d Maximum density ($nnz/(m*n)$) for $A$
 * @return Eigen::MatrixXd An $m \times n$ random binary matrix
 */
Eigen::MatrixXd randMat(unsigned int m, unsigned int n, double d) {
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, n);
  const unsigned int nnz = std::round(m * n * d);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<unsigned int> distr_row(0, m - 1);
  std::uniform_int_distribution<unsigned int> distr_col(0, n - 1);

  // we allow to draw a couple $(i, j)$ multiple times:
  // density 'd' is an upper bound for the final density of $A$
  for (unsigned int k = 0; k < nnz; ++k) {
    const unsigned int i = distr_row(gen);
    const unsigned int j = distr_col(gen);
    A(i, j) = 1;
  }

  return A;
}

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
  std::vector<unsigned int> sizes;
  constexpr double bias = 1e-5;

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

    linear.push_back(bias * n);
    quadratic.push_back(linear.back() * n);
    cubic.push_back(quadratic.back() * n);
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
  plt::title("Comparison of timings -- sparse matrices");
  plt::savefig("cx_out/matmatCOO_comparison_sparse.eps");

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

    linear.push_back(bias * n);
    quadratic.push_back(linear.back() * n);
    cubic.push_back(quadratic.back() * n);
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
  plt::title("Comparison of timings -- sparse/dense matrices");
  plt::savefig("cx_out/matmatCOO_comparison_sparsedense.eps");
}
