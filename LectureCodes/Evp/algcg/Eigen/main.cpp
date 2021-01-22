///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting algcg.m to C++/Eigen.
/// (C) 2021 SAM, D-MATH
/// Author(s): William Andersson, Vivienne Langen
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "algcg.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

using namespace Eigen;
typedef SparseMatrix<double>
    SpMat;  // declares column-major (default) sparse matrix
typedef Triplet<double> Trip;

int main() {
  // Input values.
  // Assume A is square (so evalA returs a vector size n).
  // CG requires s.p.d matrix -> use 20x20 tridiagonal matrix with diagonal=2,
  //    both of diagonals=1
  //    sparse matrix format
  // because we need cond_2(A) << 1
  unsigned int n = 20;
  unsigned int estm_entries =
      3 * n - 2;  // diagonal size n, off-diagonal size n-1
  SpMat A(n, n);

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

  VectorXd b(n);
  VectorXd x(n);
  double tol = 1e-6;
  unsigned int maxit;

  // !!  OPEN QUESTION !!
  // regarding: "..should just be operator(), no more evalA.."
  auto evalA = [A](VectorXd x) { return A * x; };

  // randome assignment of values for b,x
  // reasonable initial guess x0? or argue that parabula is bound and thus
  // converges globally
  for (unsigned int i = 0; i < n; ++i) {
    b[i] = i % 2;
    x[i] = i % 3;
  }
  // in theory the CG algorithm converges after n steps, where n=size of matrix
  // (in case for squared matrix)
  //   otherwise, but reason found in scrip?
  maxit = 22;
  // calculating and printing the solution
  VectorXd x_approx = cg(evalA, b, x, tol, maxit);
  // output
  std::cout << "The solution to A*x=b using the conjugate gradiant method is "
            << "\n";
  std::cout << " x_approx = "
            << "\n";
  std::cout << " "
            << "\n";
  std::cout << x_approx << "\n";
}
