#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <iostream>

#include "polardecomposition.hpp"

int main() {
  constexpr unsigned int m = 7;
  constexpr unsigned int n = 5;
  constexpr unsigned int k = 3;

  {
    std::cout << "Test of polar decomposition" << std::endl;
    Eigen::MatrixXd X(m, n);
    for (unsigned int i = 0; i < m; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
        X(i, j) = 1.0 / (1.0 + 2.0 * i + j);
      }
    }

    PolarDecomposition PDB(X);
    Eigen::MatrixXd Q{Eigen::MatrixXd::Identity(n, n)};
    if (PDB.Qcols() > 0) {
      PDB.applyQ(Q);
      std::cout << "Q is " << Q.rows() << " x " << Q.cols() << ":\n " << Q
                << std::endl;
    } else {
      std::cout << "Initialize Q factor before calling applyQ()!" << std::endl;
    }
    Eigen::MatrixXd M{Eigen::MatrixXd::Identity(n, n)};
    if (PDB.Mcols() > 0) {
      PDB.applyM(M);
      std::cout << "M is " << M.rows() << " x " << M.cols() << ":\n"
                << M << std::endl;
    } else {
      std::cout << "Initialize M factor before calling applyM()!" << std::endl;
    }
    if (X.rows() == Q.rows() && X.cols() == M.cols()) {
      std::cout << "|X-QM| = " << (X - Q * M).norm() << std::endl;
    } else {
      std::cout << "Shape mismatch between X and QM!" << std::endl;
    }
  }

  {
    std::cout << "\nTest of low-rank polar decomposition" << std::endl;
    Eigen::MatrixXd A(m, k);
    for (unsigned int i = 0; i < m; ++i) {
      for (unsigned int j = 0; j < k; ++j) {
        A(i, j) = 2.0 * i + j;
      }
    }
    Eigen::MatrixXd B(n, k);
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < k; ++j) {
        B(i, j) = 2.0 * i - j;
      }
    }
    PolarDecomposition PDAB(A, B);
    Eigen::MatrixXd Q{Eigen::MatrixXd::Identity(n, n)};
    if (PDAB.Qcols() > 0) {
      PDAB.applyQ(Q);
      std::cout << "Q is " << Q.rows() << " x " << Q.cols() << ":\n " << Q
                << std::endl;
    } else {
      std::cout << "Initialize Q factor before calling applyQ()!" << std::endl;
    }
    Eigen::MatrixXd M{Eigen::MatrixXd::Identity(n, n)};
    if (PDAB.Mcols() > 0) {
      PDAB.applyM(M);
      std::cout << "M is " << M.rows() << " x " << M.cols() << ":\n"
                << M << std::endl;
    } else {
      std::cout << "Initialize M factor before calling applyM()!" << std::endl;
    }
    if (A.rows() == Q.rows() && B.rows() == M.cols()) {
      std::cout << "|AB^T-QM| = " << (A * B.transpose() - Q * M).norm()
                << std::endl;
    } else {
      std::cout << "Shape mismatch between AB^T and QM!" << std::endl;
    }
  }
  return 0;
}
