///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
/* SAM_LISTING_BEGIN_0 */
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using std::cout;
using std::endl;

int main(){	
	const int n = 10;
	// Initialize a special invertible matrices
	const MatrixXd mat = MatrixXd::Identity(n,n) +
	VectorXd::Constant(n,1.0)*RowVectorXd::Constant(n,1.0);
	cout << "Matrix mat = " << endl << mat << endl;
	// Multiple right hand side vectors stored in matrix, cf. \matlab{}
	const MatrixXd B = MatrixXd::Random(n,2);
	// Solve linear system using various decompositions
	const MatrixXd X = mat.lu().solve(B);
	const MatrixXd X2 = mat.fullPivLu().solve(B);
	const MatrixXd X3 = mat.householderQr().solve(B);
	const MatrixXd X4 = mat.llt().solve(B);
	const MatrixXd X5 = mat.ldlt().solve(B);
	cout << "|X2-X| = " << (X2-X).norm() << endl;
	cout << "|X3-X| = " << (X3-X).norm() << endl;
	cout << "|X4-X| = " << (X4-X).norm() << endl;
	cout << "|X5-X| = " << (X5-X).norm() << endl;
}
/* SAM_LISTING_END_0 */
