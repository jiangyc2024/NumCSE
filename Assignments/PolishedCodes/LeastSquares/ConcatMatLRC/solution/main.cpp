/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "concatmatlrc.hpp"

int main() {
  std::cout << "Economical SVD of concatenated matrices" << std::endl;
  {
    // Test data
    constexpr unsigned int m = 7;
    constexpr unsigned int n = 5;
    constexpr unsigned int k = 2;
    std::cout << "CASE I" << std::endl;
    std::cout << "#################################" << std::endl;
    std::cout << "Printing the matrices A1,B1,A2,B2 \n" << std::endl;
    Eigen::MatrixXd A1 =
        (Eigen::MatrixXd(m, k) << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
            .finished();
    Eigen::MatrixXd B1 =
        (Eigen::MatrixXd(n, k) << 1, -2, 3, -4, 5, -6, 7, -8, 9, 10).finished();
    Eigen::MatrixXd A2 = (Eigen::MatrixXd(m, k) << -1, 2, 3, -4, 5, 6, 7, -8, 9,
                          10, 11, 12, 13, 14)
                             .finished();
    Eigen::MatrixXd B2 =
        (Eigen::MatrixXd(n, k) << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10).finished();
    std::cout << "A1 = " << std::endl << A1 << std::endl;
    std::cout << "\nB1 = " << std::endl << B1 << std::endl;
    std::cout << "\nA2 = " << std::endl << A2 << std::endl;
    std::cout << "\nB2 = " << std::endl << B2 << std::endl;
    auto [U, s, V] = eco_svd_concat(A1, B1, A2, B2);
    // std::cout << "\nSingular values = " << s.transpose() << std::endl;
    // std::cout << "\nFactor U = " << std::endl << U << std::endl;
    // std::cout << "\nFactor V = " << std::endl << V << std::endl;
    // std::cout << "\ntest_eco_svd_concat() result: "
    //           << (test_eco_svd_concat(A1, B1, A2, B2) ? "success\n"
    //                                                   : "failure\n")
    //           << std::endl;
  }

  {
    constexpr unsigned int m = 200;
    constexpr unsigned int n = m / 2;
    constexpr unsigned int k = 20;
    Eigen::MatrixXd A1(m, k);
    Eigen::MatrixXd A2(m, k);
    Eigen::MatrixXd B1(n, k);
    Eigen::MatrixXd B2(n, k);

    for (unsigned int i = 0; i < m; ++i) {
      for (unsigned int j = 0; j < k; ++j) {
        A1(i, j) = 1.0 / (i + j + 1.0);
        A2(i, j) = 1.0 / (i + j + 3.0);
      }
    }
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < k; ++j) {
        B1(i, j) = 1.0 / (i + j + 2.0);
        B2(i, j) = 1.0 / (i + j + 4.0);
      }
    }
    std::cout << "CASE II" << std::endl;
    std::cout << "#################################" << std::endl;
    std::cout << "test_eco_svd_concat() result: "
              << (test_eco_svd_concat(A1, B1, A2, B2) ? "success\n"
                                                      : "failure\n")
              << std::endl;

    std::cout << "Low rank approximation result \n" << std::endl;
    double tol = 0.5;
    for (unsigned int k = 0; k < 6; ++k) {
      auto [AM, BM] = concat_low_rank_best(A1, B1, A2, B2, tol);
      Eigen::MatrixXd X(m, 2 * n);
      X.block(0, 0, m, n) = A1 * B1.transpose();
      X.block(0, n, m, n) = A2 * B2.transpose();
      std::cout << "tol = " << tol << ": rank = " << AM.cols()
                << ", approximation error = "
                << (AM * BM.transpose() - X).norm() << std::endl;
      tol /= 10.0;
    }
  }
  return 0;
}
