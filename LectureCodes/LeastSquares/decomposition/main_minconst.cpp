///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <iostream>

#include "./minconst.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

//NOLINTBEGIN (bugprone-exception-escape)
int main() {
  // Initialize a random matrix
  const Eigen::Index m = 9;
  const Eigen::Index n = 4;
  MatrixXd A(m,n);
  for(int i=0;i<m;i++) {
    for(int j=0;j<n;j++) {
      A(i,j) = static_cast<double>(i-j)/(i+j+1);
    }
  }

  cout << "rank(A) = " << A.jacobiSvd().rank() << endl;
  cout << "Minimizer of x -> |Ax| on unit sphere:" << endl;
  VectorXd x(n); 
  const double smin = minconst::minconst(x,A);
  cout << "x = " << x.transpose() << ", minimum = " << smin << endl;
  return 0;
}
//NOLINTEND
  
