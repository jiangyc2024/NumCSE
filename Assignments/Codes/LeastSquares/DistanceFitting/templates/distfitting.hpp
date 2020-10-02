#pragma once

#include "totriplets.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace Eigen;


SparseMatrix<double> initA(unsigned int n) {
  unsigned int rows = n * (n - 1) / 2;
  unsigned int cols = n - 1;
  SparseMatrix<double> A(rows, cols);
  // TO DO (0-1.b)
  // START
  
  // END
  A.makeCompressed();
  return A;
}



VectorXd solveExtendedNormalEquations(const MatrixXd &D) {
  VectorXd x;
  unsigned int n = D.cols();
  unsigned int m = n * (n - 1) / 2;
  // TO DO (0-1.c)
  // START
  
  // END
  return x;
}



VectorXd solveNormalEquations(const MatrixXd &D) {
  VectorXd x;
  // TO DO (0-1.e)
  // START
  
  // END
  return x;
}
