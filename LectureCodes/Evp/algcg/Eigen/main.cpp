///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting algcg.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
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
typedef SparseMatrix<double> SpMat; // declares column-major (default) sparse matrix   
typedef Triplet<double> Trip; 

int main() {
  // Input values.
  // Assume A is square (so evalA returs a vector size n).
  // CG requires s.p.d matrix -> use 20x20 tridiagonal matrix with diagonal=2, 
  //    both of diagonals=1 
  //    sparse matrix format 
  unsigned int n = 20;
  unsigned int estm_entries = 3*n-2; // diagonal size n, off-diagonal size n-1 
  SpMat A(n,n);

  // randome access of sparse matrix is expensive 
  // -> first build list of triplets, then convert to sparse matrix 
  std::vector<Trip> tripletList;
  tripletList.reserve(estm_entries); 
  for(unsigned int i=0; i<n; ++i){
    for(unsigned int j=0; j<n; ++j){
        // diagonals
        if(i==j){ tripletList.push_back(Trip(i,j,2));}
        // off-diagonals
        if((i-j)==1){ tripletList.push_back(Trip(i,j,1));}  //when internet, figure out why abs(i-j)==1 isn't working
        if((j-i)==1){ tripletList.push_back(Trip(i,j,1));}  //when internet, figure out why abs(i-j)==1 isn't working
    }
  }
  A.setFromTriplets(tripletList.begin(), tripletList.end());

  std::cout << A;

  VectorXd b(n);
  VectorXd x(n);
  double tol = 1e-6;
  unsigned int maxit;


  auto evalA = [A](VectorXd x) { return A * x; };

  // let's try randome values for b,x 
  for(unsigned int i=0; i<n; ++i){
    b[i] = i;
    x[i] = n-i;
  }

  maxit = 5;

  VectorXd x_approx = cg(evalA, b, x, tol, maxit);

  std::cout << x_approx << "\n";
}
