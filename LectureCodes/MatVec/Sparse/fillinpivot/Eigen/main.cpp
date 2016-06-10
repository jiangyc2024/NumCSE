#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <figure/figure.hpp>
using namespace std;
using namespace Eigen;

int main () {
	#pragma begin<0>
	// Study of fill-in with LU-factorization due to pivoting
	MatrixXd A(11,11); A.setZero();
	A.diagonal() = VectorXd::LinSpaced(11,1,11).cwiseInverse();
	A.col(10).setConstant(2); A.row(10).setConstant(2);
	auto solver = A.lu();
	MatrixXd L = MatrixXd::Identity(11,11);
	L += solver.matrixLU().triangularView<StrictlyLower>();
	MatrixXd U = solver.matrixLU().triangularView<Upper>();
	// Plotting
	mgl::Figure fig1, fig2, fig3, fig4;
	fig1.spy(A);	fig1.setFontSize(4); fig1.grid(false);
	fig1.title("arrow matrix A");	fig1.save("fillinpivotA");
	fig2.spy(L);	fig2.setFontSize(4); fig2.grid(false);
	fig2.title("L factor");	fig2.save("fillinpivotL");
	fig3.spy(U);	fig3.setFontSize(4); fig3.grid(false);
	fig3.title("U factor");	fig3.save("fillinpivotU");
	std::cout  << A << std::endl;
	#pragma end<0>
	return 0;
}
