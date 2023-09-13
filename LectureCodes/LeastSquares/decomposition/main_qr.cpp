///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include "./decomp.hpp"
#include <Eigen/Dense>
#include <iostream>


using Eigen::MatrixXd;
using std::cout;
using std::endl;

int main() {
  // Initialize a random matrix
  const MatrixXd A = MatrixXd::Random(9,4);
  // Matrices storing QR-factors
  MatrixXd R; // R-factors
  MatrixXd R0; 
  MatrixXd Q; // (thin) Q-factor
  MatrixXd Q0; 

  cout << "full QR-decompositoion of A = " << endl << A << endl;
  std::tie(Q,R) = decomp::qr_decomp_full(A);
  cout << "Q of size " << Q.rows() << "x" << Q.cols() << " = " << endl << Q << endl;
  cout << "R of size " << R.rows() << "x" << R.cols() << " = " << endl << R << endl;
  cout << "|A-QR| = " << (A-Q*R).norm() << endl << endl;

  cout << "economical QR-decompositoion of A = " << endl << A << endl;
  std::tie(Q0,R0) = decomp::qr_decomp_eco(A);
  cout << "Q0 of size " << Q0.rows() << "x" << Q0.cols() << " = " << endl << Q0 << endl;
  cout << "R0 of size " << R0.rows() << "x" << R0.cols() << " = " << endl << R0 << endl;
  cout << "|A-QR| = " << (A-Q0*R0).norm() << endl;

  return 0;
}
