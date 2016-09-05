#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

int main(){	
	int n = 10;
	// Initialize a special invertible matrices
	MatrixXd mat = MatrixXd::Identity(n,n) +
	 VectorXd::Constant(n,1.0)*RowVectorXd::Constant(n,1.0);
	cout << "Matrix mat = " << endl << mat << endl;
	// Multiple right hand side vectors stored in matrix, cf. \matlab{}
	MatrixXd B = MatrixXd::Random(n,2);
	// Solve linear system using various decompositions
	MatrixXd X = mat.lu().solve(B);
	MatrixXd X2 = mat.fullPivLu().solve(B);
	MatrixXd X3 = mat.householderQr().solve(B);
	MatrixXd X4 = mat.llt().solve(B);
	MatrixXd X5 = mat.ldlt().solve(B);
	cout << "|X2-X| = " << (X2-X).norm() << endl;
	cout << "|X3-X| = " << (X3-X).norm() << endl;
	cout << "|X4-X| = " << (X4-X).norm() << endl;
	cout << "|X5-X| = " << (X5-X).norm() << endl;
}
