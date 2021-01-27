// this is the c++ equivalent code of example LSE for Ex. 10.3.0.11 on page 630
// MATLAB translated code 

///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting pcgbase.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): Vivienne Langen
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "pcgbase.hpp"

#include <iostream>
#include <Eigen/Dense> 
#include <Eigen/Sparse>
#include <Eigen/Core>

using namespace Eigen;
typedef SparseMatrix<double> SpMat;  // declares column-major (default) sparse matrix
typedef Triplet<double> Trip;
typedef std::tuple<VectorXd, VectorXd> tuple2;

/* Example 10.3.0.11 for LSE with tridiagonal preconditioning */
int main() { 
    // input values 
    // matrix size
    unsigned int n = 20; // has to be even number if example of script is used 
    double tol = 1e-4;
    unsigned int maxit = 1000;
    // pcgbase method requires s.p.d matrix, suitable for preconditioning   
    unsigned int estm_entries = n + 2 * (n - 1) + 2*(n/2) ;  // diagonal size n, 2 times off-diagonal size n-1, 2 times diagonal at n/2
    SpMat A(n, n);
    VectorXd b(n); // r.h.s.
    VectorXd x0(n); // initial guess
    // randome access of sparse matrix is expensive
    // -> first build list of triplets, then convert to sparse matrix
    std::vector<Trip> tripletList;
    tripletList.reserve(estm_entries);
    // filling matrix with diagonals=2, off-diagonals=1
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
        // diagonals
        if (i == j) {
          tripletList.push_back(Trip(i, j, 2+2/float(n) ));
          }
        // off-diagonals
        // convert to 'int' to use abs()
        int dd = i - j;
        if (abs(dd) == 1) {
           tripletList.push_back(Trip(i, j, -1));
           }
        if (abs(dd) == n/2) {
           tripletList.push_back(Trip(i, j, 1/float(n) ));
           }
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    // according to example
    b = x0 = VectorXd::Ones(n);
    // !!  OPEN QUESTION !!
    // regarding: "..should just be operator(), no more evalA.."
    auto evalA = [A](VectorXd x) { return A * x; };
    auto invB = [](VectorXd x) { return x; };
    // solution vectors
    VectorXd r; //residual
    VectorXd x; //approximate solution 
    tuple2 sol = std::make_tuple(x, r);
    sol  = pcgbase<tuple2> (evalA, invB, b, x0, tol, maxit);
}
