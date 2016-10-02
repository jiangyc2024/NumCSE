///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <Eigen/Dense>

#include "./svdspaces.hpp"

using namespace Eigen;
using namespace std;

int main(void) {
  // Initialize a random matrix
  const int m_init=19,n_init=14,r=6;

  MatrixXd A = MatrixXd::Zero(m_init,n_init);
  for (int j =0; j<r; j++) A += VectorXd::Random(m_init)*RowVectorXd::Random(n_init);

  cout << "rank_eigen(A) = " << rank_eigen(A) << endl;
  cout << "rank_ncse(A) = " << rank_ncse(A) << endl;
  MatrixXd Z = nullspace(A);
  int dz = Z.cols();
  cout << "ONB of kernel of A (dim = " << dz << ") = " << endl << Z << endl;
  cout << "|I-Z^T*Z| = " << (MatrixXd::Identity(dz,dz)-Z.transpose()*Z).norm() << endl;
  cout << "|A*Z| = " << (A*Z).norm() << endl;

  MatrixXd B = rangespace(A);
  int r = B.cols();
  cout << "ONB of range space of A (dim = " << r << ") = " << endl << B << endl;
  
  return 0;
}
