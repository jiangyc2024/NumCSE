/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "auxiliary_functions.hpp"
#include "gaussseidelcrs.hpp"

int main() {
  std::cout << "C++ program for course Numerical Methods for CSE" << std::endl;
  std::cout << "Gauss-Seidel iteration with CRS matrix \n" << std::endl;

  CRSMatrix M;
  M.m = 6;
  M.n = 6;
  M.val = {10.0, -2.0, 3.0, 9.0, 3.0,  2.0, 8.0, 2.0, 3.0, 1.0,
           10.0, 5.0,  2.0, 3.0, 13.0, 2.0, 4.0, 2.0, 11.0};
  M.col_ind = {0, 4, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5};
  M.row_ptr = {0, 2, 5, 8, 12, 16, 19};

  Eigen::MatrixXd A(6, 6);
  // clang-format off
  A << 10, 0, 0, 0, -2, 0,
       3, 9, 0, 0, 0, 3,
       0, 2, 8, 2, 0, 0,
       3, 0, 1, 10, 5, 0,
       0, 2, 0, 3, 13, 2,
       0, 4, 0, 0, 2, 11;
  // clang-format on

  {
    std::cout << "Gauss-Seidel step for dense matrix (available for reference)"
              << std::endl;
    Eigen::VectorXd x{Eigen::VectorXd::LinSpaced(6, 1.0, 6.0)};
    Eigen::VectorXd b{Eigen::VectorXd::LinSpaced(6, 1.0, 6.0)};
    std::cout << (GaussSeidelstep_generic(A, b, x) ? "Success" : "Failure")
              << std::endl;
    std::cout << "x = " << x.transpose() << std::endl;
  }
  {
    std::cout << "\nGauss-Seidel step for CRS matrix (should match the dense "
                 "version above)"
              << std::endl;
    Eigen::VectorXd x{Eigen::VectorXd::LinSpaced(6, 1.0, 6.0)};
    Eigen::VectorXd b{Eigen::VectorXd::LinSpaced(6, 1.0, 6.0)};
    std::cout << (GaussSeidelstep_crs(M, b, x) ? "Success" : "Failure")
              << std::endl;
    std::cout << "x = " << x.transpose() << std::endl;
  }
  {
    std::cout << "\nGauss-Seidel iteration for CRS matrix (should return "
                 "vector of zeros)"
              << std::endl;
    Eigen::VectorXd x{Eigen::VectorXd::LinSpaced(6, 1.0, 6.0)};
    Eigen::VectorXd b{Eigen::VectorXd::Zero(6)};
    std::cout << (GaussSeidel_iteration(M, b, x) ? "Success" : "Failure")
              << std::endl;
    std::cout << "x = " << x.transpose() << std::endl;
  }
  return 0;
}
