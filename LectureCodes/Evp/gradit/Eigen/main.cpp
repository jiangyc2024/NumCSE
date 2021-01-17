///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting gradit.m to C++/Eigen.
/// (C) 2021 SAM, D-MATH
/// Author(s): William Andersson, Vivienne Langen
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "gradit.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

using namespace Eigen;
typedef SparseMatrix<double>
    SpMat;  // declares column-major (default) sparse matrix
typedef Triplet<double> Trip;

int main() {
  // Input values A, b, x
  // Gradient method requires s.p.d, symmetric matrix -> use 20x20 tridiagonal
  // matrix
  //    with diagonal=2, both of diagonals=1
  //    sparse matrix format
  unsigned int n = 20;
  // estimated entries for matrix
  unsigned int estm_entries =
      3 * n - 2;  // diagonal size n, 2 times off-diagonal size n-1
  SpMat A(n, n);
  VectorXd b(n);
  VectorXd x0(n);
  // randome access of sparse matrix is expensive
  // -> first build list of triplets, then convert to sparse matrix
  std::vector<Trip> tripletList;
  tripletList.reserve(estm_entries);
  // filling matrix with diagonals=2, off-diagonals=1
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      // diagonals
      if (i == j) {
        tripletList.push_back(Trip(i, j, 2));
      }
      // off-diagonals
      // convert to 'int' to use abs()
      int dd = i - j;
      if (abs(dd) == 1) {
        tripletList.push_back(Trip(i, j, 1));
      }
    }
  }
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  // where do these come from??
  double rtol = 1e-6;
  double atol = 1e-6;
  // find reasonable maxit?
  unsigned int maxit = 20;
  // !!  OPEN QUESTION !!
  // regarding: "..should just be operator(), no more evalA.."
  auto evalA = [A](VectorXd x) { return A * x; };
  // randome assignment of values for b,x
  // reasonable initial guess x0? or argue that parabula is bound and thus
  // converges globally
  for (unsigned int i = 0; i < n; ++i) {
    b[i] = i % 2;
    x0[i] = i % 3;
  }
  // Output
  VectorXd x_approx = gradit(evalA, A, b, x0, rtol, atol, maxit);
  std::cout << "The solution to A*x=b using the gradiant method with steepest "
               "descent is "
            << "\n";
  std::cout << " x_approx = "
            << "\n";
  std::cout << " "
            << "\n";
  std::cout << x_approx << "\n";
  std::cout << x_approx << "\n";
}
