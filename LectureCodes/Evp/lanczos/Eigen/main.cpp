#include <iostream>
#include <Eigen/Dense>
#include <figure.hpp>
#include "lanczos.hpp"
#include "lanczosev.hpp"

#include <Eigen/Eigenvalues>
#include <eigensolversort.hpp>

int main()
{
	// create test matrix
	
	int n = 100;
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n,n);
	M.diagonal(-1) = -0.5*Eigen::VectorXd::Ones(n-1);
	M.diagonal(0) =   2.*Eigen::VectorXd::Ones(n);
	M.diagonal(1) = -1.5*Eigen::VectorXd::Ones(n-1);

	auto qr = M.householderQr();   // generate orthogonal matrix
	Eigen::MatrixXd Q = qr.householderQ();

	// eigenvalues \Blue{$1,2,3,\ldots,100$}
	Eigen::VectorXd diag= Eigen::VectorXd::LinSpaced(n,1,n);
	Eigen::MatrixXd A = Q.transpose()*diag.asDiagonal()*Q;

	// test lancos
	int k = 30;
	Eigen::VectorXd z0 = Eigen::VectorXd::Ones(n);
	Eigen::VectorXd alpha;
	Eigen::VectorXd beta;
	Eigen::MatrixXd V;
	lanczos(A, k, z0, V, alpha, beta);

	std::cout << "alpha:" << std::endl;
	std::cout << alpha << std::endl;
	std::cout << "beta:" << std::endl;
	std::cout << beta << std::endl;


	// test lancosev 
	Eigen::EigenSolver<Eigen::MatrixXd> esolver(A);
	Eigen::MatrixXd U; 
	Eigen::VectorXd ev; 
	std::tie(ev, U) = eigensolversort(esolver);

	
	int j = 5;
	Eigen::VectorXd ones = Eigen::VectorXd::Ones(100);
	Eigen::MatrixXd res = lanczosev(A,k,j,ones);

	{
		mgl::Figure fig;
		Eigen::VectorXd y1 = res.col(j) - Eigen::MatrixXd::Constant(k, 1, ev(n-1));
		Eigen::VectorXd y2 = res.col(j-1) - Eigen::MatrixXd::Constant(k, 1, ev(n-2));
		Eigen::VectorXd y3 = res.col(j-2) - Eigen::MatrixXd::Constant(k, 1, ev(n-3));
		Eigen::VectorXd y4 = res.col(j-3) - Eigen::MatrixXd::Constant(k, 1, ev(n-4));
		Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(k,1,k);
		fig.setlog(false, true);
		fig.plot(x, 		  y1.cwiseAbs(), 		   "r+").label("lamda n");
		fig.plot(x.tail(k-1), y2.cwiseAbs().tail(k-1), "m-^").label("lamda n-1");
		fig.plot(x.tail(k-2), y3.cwiseAbs().tail(k-2), "b-o").label("lamda n-2");
		fig.plot(x.tail(k-3), y4.cwiseAbs().tail(k-3), "k-*").label("lamda n-3");
		
		fig.legend(1,1);
		fig.xlabel("step of Lanzcos process");
		fig.ylabel("|Ritz value-eigenvalue|");
		fig.setFontSize(3);
		fig.save("lanczosev1");
	}
}
