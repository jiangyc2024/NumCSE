///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <iostream>
# include <Eigen/Dense>
# include "./vandermonde.hpp"

using namespace Eigen;
using namespace std;

int main() {
  const int n = 7;
  VectorXd t = VectorXd::LinSpaced(n,1,n);
  cout << "t = " << t.transpose() << endl << "V = " << endl << vander(t) << endl;
  return 0;
}
