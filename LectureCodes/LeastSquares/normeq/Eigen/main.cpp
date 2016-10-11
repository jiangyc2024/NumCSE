///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <limits>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "normeeqsolve.hpp"

int main(void) {
  // Dimensions of sample least squares problem
  int m=7,n=3; MatrixXd A(m,n); VectorXd b(m);
  // Initialization of testing matrix
  for(int j=0;j<m;j++) {
    b(j) = 1./(j+1);  for(int i=0;i<n;i++) A(j,i) = 1+2*i*i-j;
  }
  A(0,0) = 0;
  // Output of overdetermined system for testing
  cout << "system matrix A = " << endl << A << endl << "r.h.s. b = " << b << endl;
  // Solve normal equations
  VectorXd x = normeqsolve(A,b);
  // Output presumable least squares solution
  cout << "x = " << x.transpose() << endl;
  return 0;
}
