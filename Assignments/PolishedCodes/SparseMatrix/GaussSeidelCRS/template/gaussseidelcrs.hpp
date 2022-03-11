/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */
#ifndef GAUSSSEIDELCRSHPP
#define GAUSSSEIDELCRSHPP

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

#define _USE_MATH_DEFINES

struct CRSMatrix {
  unsigned int m;                     // number of rows
  unsigned int n;                     // number of columns
  std::vector<double> val;            // value array
  std::vector<unsigned int> col_ind;  // same length as value array
  std::vector<unsigned int> row_ptr;  // length m+1, row_ptr[m] == val.size()
};

// One step of a Gauss-Seidel iteration for the linear system of equations Ax=b
// with A given in a raw CRS format
bool GaussSeidelstep_crs(const CRSMatrix &A, const Eigen::VectorXd &b,
                         Eigen::VectorXd &x) {
  assert(A.n == A.m && "Matrix must be square");
  assert(A.n == b.size() && "Vector b length mismatch");
  assert(A.n == x.size() && "Vector x length mismatch");

  bool out;
  // TODO : Write a function that executes one step of the Gauss-Seidel 
  // iteration. It should return false if the iteration is not well defined, 
  // otherwise it returns true.
  // START 
  
  // END 
  return out;
}

// Gauss-Seidel iteration for matrix A and vector b with correction-based
// termination
bool GaussSeidel_iteration(const CRSMatrix &A, const Eigen::VectorXd &b,
                           Eigen::VectorXd &x, double atol = 1.0E-8,
                           double rtol = 1.0E-6, unsigned int maxit = 100) {
  assert(A.n == A.m && "Matrix must be square");
  assert(A.n == b.size() && "Vector b length mismatch");
  assert(A.n == x.size() && "Vector x length mismatch");

  bool out;
  // TODO : Write a function that carries out the Gauss-Seidel iteration for
  // the given inputs A,b and initial guess x^0 passed in x. x should also 
  // return the final iterate. If the iteration doesn't converge in maxit steps
  // the function should return false, otherwise it returns true.
  // START
  
  // END

  return out;
}

#endif
