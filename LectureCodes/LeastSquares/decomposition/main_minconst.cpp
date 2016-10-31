///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <Eigen/Dense>

#include "./minconst.hpp"

using namespace Eigen;
using namespace std;

int main(void) {
  // Initialize a random matrix
  MatrixXd::Index m=9,n=4;
  MatrixXd A(m,n);
  for(int i=0;i<m;i++) for(int j=0;j<n;j++) A(i,j) = ((double)i-j)/(i+j+1);

  cout << "rank(A) = " << A.jacobiSvd().rank() << endl;
  cout << "Minimizer of x -> |Ax| on unit sphere:" << endl;
  VectorXd x(n); double smin = minconst(x,A);
  cout << "x = " << x.transpose() << ", minimum = " << smin << endl;
  return 0;
}
  
