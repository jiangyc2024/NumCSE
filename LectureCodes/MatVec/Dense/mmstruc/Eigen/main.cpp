#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

int main () {
	int n = 100;
	MatrixXd A(n,n), B(n,n); A.setZero(); B.setZero();
	A.diagonal() = VectorXd::LinSpaced(n,1,n);
	A.col(n-1) = VectorXd::LinSpaced(n,1,n);
	A.row(n-1) = RowVectorXd::LinSpaced(n,1,n);
	B = A.colwise().reverse();
	MatrixXd C = A*A, D = A*B;
	mgl::Figure fig1, fig2, fig3, fig4;
	fig1.spy(A);	fig1.save("Aspy_cpp");
	fig2.spy(B);	fig2.save("Bspy_cpp");
	fig3.spy(C);	fig3.save("Cspy_cpp");
	fig4.spy(D);	fig4.save("Dspy_cpp");
	return 0;
}
