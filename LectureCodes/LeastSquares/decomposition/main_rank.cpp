///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <iostream>

#include "./svdspaces.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using std::cout;
using std::endl;

int main() {
    // Initialize a random matrix
    const int m_init = 19;
    const int n_init = 14;
    const int r = 6;

    MatrixXd A = MatrixXd::Zero(m_init,n_init);
    for (int j =0; j<r; j++) {
        A += VectorXd::Random(m_init)*RowVectorXd::Random(n_init);
    }
    cout << "rank_eigen(A) = " << rank_eigen(A) << endl;
    cout << "rank_ncse(A) = " << rank_by_svd(A) << endl;
    MatrixXd Z = nullspace(A);
    const Eigen::Index dz = Z.cols();

    cout << "ONB of kernel of A (dim = " << dz << ") = " << endl << Z << endl;
    cout << "|I-Z^T*Z| = " << (MatrixXd::Identity(dz,dz)-Z.transpose()*Z).norm() << endl;
    cout << "|A*Z| = " << (A*Z).norm() << endl;

    const MatrixXd B = rangespace(A);
    cout << "ONB of range space of A (dim = " << r << ") = " << endl << B.cols() << endl;

    return 0;
}
