///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <iostream>

#include "./decomp.hpp"
#include "./svdspaces.hpp"

using Eigen::MatrixXd;
using std::cout;
using std::endl;
using std::min;

int main() {
  // Initialize a random matrix
  const int m_init = 9;
  const int n_init = 4;
  MatrixXd A = MatrixXd::Random(m_init, n_init);
  
  // SVD factors
  MatrixXd U; 
  MatrixXd Sigma;
  MatrixXd V; 

  {
    cout << "full SVD of A = " << endl << A << endl;
    const Eigen::Index m = A.rows();
    const Eigen::Index n = A.cols();
    std::tie(U, Sigma, V) = decomp::svd_full(A);
    cout << "U of size " << U.rows() << "x" << U.cols() << " = " << endl
         << U << endl;
    cout << "S of size " << Sigma.rows() << "x" << Sigma.cols() << " = " << endl
         << Sigma << endl;
    cout << "V of size " << V.rows() << "x" << V.cols() << " = " << endl
         << V << endl;
    cout << "|I-U^T*U|" << (MatrixXd::Identity(m, m) - U.transpose() * U).norm()
         << endl;
    cout << "|I-V^T*U|" << (MatrixXd::Identity(n, n) - V.transpose() * V).norm()
         << endl;
    cout << "|A-U*S*V| = " << (A - U * Sigma * V.transpose()).norm() << endl
         << endl;
  }

  {
    cout << "Economical SVD of A = " << endl << A << endl;
    const Eigen::Index m = A.rows();
    const Eigen::Index n = A.cols();
    const Eigen::Index p = min(m, n);
    std::tie(U, Sigma, V) = decomp::svd_eco(A);
    cout << "U of size " << U.rows() << "x" << U.cols() << " = " << endl
         << U << endl;
    cout << "S of size " << Sigma.rows() << "x" << Sigma.cols() << " = " << endl
         << Sigma << endl;
    cout << "V of size " << V.rows() << "x" << V.cols() << " = " << endl
         << V << endl;
    cout << "|I-U^T*U|" << (MatrixXd::Identity(p, p) - U.transpose() * U).norm()
         << endl;
    cout << "|I-V^T*U|" << (MatrixXd::Identity(p, p) - V.transpose() * V).norm()
         << endl;
    cout << "|A-U*S*V| = " << (A - U * Sigma * V.transpose()).norm() << endl
         << endl;
  }

  const MatrixXd AT = A.transpose();
  {
    cout << "Economical SVD of A^T = " << endl << AT << endl;
    const Eigen::Index m = AT.rows();
    const Eigen::Index n = AT.cols();
    const Eigen::Index p = min(m, n);
    std::tie(U, Sigma, V) = decomp::svd_eco(AT);
    cout << "U of size " << U.rows() << "x" << U.cols() << " = " << endl
         << U << endl;
    cout << "S of size " << Sigma.rows() << "x" << Sigma.cols() << " = " << endl
         << Sigma << endl;
    cout << "V of size " << V.rows() << "x" << V.cols() << " = " << endl
         << V << endl;
    cout << "|I-U^T*U|" << (MatrixXd::Identity(p, p) - U.transpose() * U).norm()
         << endl;
    cout << "|I-V^T*U|" << (MatrixXd::Identity(p, p) - V.transpose() * V).norm()
         << endl;
    cout << "|A-U*S*V| = " << (AT - U * Sigma * V.transpose()).norm() << endl
         << endl;
  }

  return 0;
}
