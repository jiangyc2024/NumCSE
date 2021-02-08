///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting pcgbase.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson, Vivienne Langen
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

/* Standard example */
int main() { 
    // input values 
    // matrix size
    unsigned int n = 20; // has to be even number if example of script is used 
    double tol = 1e-4;
    int maxit = 1000;
    // pcgbase method requires s.p.d matrix, suitable for preconditioning   
    unsigned int estm_entriesA = n + 2 * (n - 1) ;  // diagonal size n, 2 times off-diagonal size n-1
    unsigned int estm_entriesB = n ;  // diagonal size n
    SpMat A(n, n);
    SpMat B(n, n);
    VectorXd b(n); // r.h.s.
    VectorXd x0(n); // initial guess
    // randome access of sparse matrix is expensive
    // -> first build list of triplets, then convert to sparse matrix
    std::vector<Trip> tripletListA;
    std::vector<Trip> tripletListB;
    tripletListA.reserve(estm_entriesA);
    tripletListB.reserve(estm_entriesB);
    // filling matrix A with diagonals=1,2,3,4,..,n ,off-diagonals=1
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
        // diagonals
        if (i == j) {
          tripletListA.push_back(Trip(i, j, i+1));
        }
        // off-diagonals
        // convert to 'int' to use abs()
        int dd = i - j;
        if (abs(dd) == 1) {
           tripletListA.push_back(Trip(i, j, 1));
        }
      }
    } 
    // matrix B is the inverse of the diagonal matrix with diag=1,2,3,4,..,n
    for (unsigned int i = 0; i < n; ++i) {
      // diagonals
      tripletListB.push_back(Trip(i, i, 1/float(i+1)));
    }
    A.setFromTriplets(tripletListA.begin(), tripletListA.end());
    B.setFromTriplets(tripletListB.begin(), tripletListB.end());
    // according to script example 
    b = x0 = VectorXd::Ones(n);
    // !!  OPEN QUESTION !!
    // regarding: "..should just be operator(), no more evalA.."
    auto evalA = [A](VectorXd x) { return A * x; };
    auto idB = [](VectorXd x) { return x; };
    // matrix for tridiagonal preconditioning
    auto tridiagB = [B](VectorXd x) { return B * x; };
    // solution vectors
    VectorXd r; //residual
    VectorXd x; //approximate solution 
    std::pair<VectorXd, VectorXd> sol_nop, sol_tridiagp;
    sol_nop = sol_tridiagp = std::make_pair(x, r) ;
    // call CG method with no pre-conditioning 
    sol_nop  = pcgbase(evalA, idB, b, x0, tol, maxit);
    // call CG method with tridiagnal pre-conditioning
    sol_tridiagp = pcgbase(evalA, tridiagB, b, x0, tol, maxit);
    std::cout << "Approximated solution with no pre-conditioning x = " << "\n" << std::get<0>(sol_nop) << "\n";
    std::cout << "Approximated solution with tridiagonal pre-conditioning x = " << "\n" << std::get<0>(sol_tridiagp) << "\n";
}
