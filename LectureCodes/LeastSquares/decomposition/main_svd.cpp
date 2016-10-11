///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <Eigen/Dense>

#include "./decomp.hpp"
#include "./svdspaces.hpp"

using namespace Eigen;
using namespace std;

int main(void) {
  // Initialize a random matrix
  const int m_init=9,n_init=4; 
  MatrixXd A = MatrixXd::Random(m_init,n_init);
  MatrixXd U,Sigma,V; // SVD factors

  {
    cout << "full SVD of A = " << endl << A << endl;
    const int m=A.rows(),n=A.cols();
    std::tie(U,Sigma,V) = svd_full(A);
    cout << "U of size " << U.rows() << "x" << U.cols() << " = " << endl << U << endl;
    cout << "S of size " << Sigma.rows() << "x" << Sigma.cols() << " = " << endl << Sigma << endl;
    cout << "V of size " << V.rows() << "x" << V.cols() << " = " << endl << V << endl;
    cout << "|I-U^T*U|" << (MatrixXd::Identity(m,m)-U.transpose()*U).norm() << endl;
    cout << "|I-V^T*U|" << (MatrixXd::Identity(n,n)-V.transpose()*V).norm() << endl;
    cout << "|A-U*S*V| = " << (A-U*Sigma*V.transpose()).norm() << endl << endl;
  }

  {
    cout << "Economical SVD of A = " << endl << A << endl;
    const int m=A.rows(),n=A.cols();
    const int p = min(m,n);
    std::tie(U,Sigma,V) = svd_eco(A);
    cout << "U of size " << U.rows() << "x" << U.cols() << " = " << endl << U << endl;
    cout << "S of size " << Sigma.rows() << "x" << Sigma.cols() << " = " << endl << Sigma << endl;
    cout << "V of size " << V.rows() << "x" << V.cols() << " = " << endl << V << endl;
    cout << "|I-U^T*U|" << (MatrixXd::Identity(p,p)-U.transpose()*U).norm() << endl;
    cout << "|I-V^T*U|" << (MatrixXd::Identity(p,p)-V.transpose()*V).norm() << endl;
    cout << "|A-U*S*V| = " << (A-U*Sigma*V.transpose()).norm() << endl << endl;
  }

  MatrixXd AT = A.transpose();
  {
    cout << "Economical SVD of A^T = " << endl << AT << endl;
    const int m=AT.rows(),n=AT.cols();
    const int p = min(m,n);
    std::tie(U,Sigma,V) = svd_eco(AT);
    cout << "U of size " << U.rows() << "x" << U.cols() << " = " << endl << U << endl;
    cout << "S of size " << Sigma.rows() << "x" << Sigma.cols() << " = " << endl << Sigma << endl;
    cout << "V of size " << V.rows() << "x" << V.cols() << " = " << endl << V << endl;
    cout << "|I-U^T*U|" << (MatrixXd::Identity(p,p)-U.transpose()*U).norm() << endl;
    cout << "|I-V^T*U|" << (MatrixXd::Identity(p,p)-V.transpose()*V).norm() << endl;
    cout << "|A-U*S*V| = " << (AT-U*Sigma*V.transpose()).norm() << endl << endl;
  }
    
  return 0;
}
