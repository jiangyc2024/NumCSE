#include <iostream>
#include <Eigen/Dense>
#include <figure.hpp>
#include "eigensolversort.hpp"
#include "kryleig.hpp"
#include "spdiags.hpp"


int main()
{
	// Generate matrix
	
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

	Eigen::EigenSolver<Eigen::MatrixXd> esolver(A);

	Eigen::VectorXd d0; // eigenvalues
	Eigen::MatrixXd V0; // eigenvectors
	std::tie(d0, V0) = eigensolversort(esolver);

	Eigen::VectorXd ev0(6);
	ev0.head(3) = d0.head(3);
	ev0.tail(3) = d0.tail(3);

	int N = 30-5;
	Eigen::MatrixXd errev(N,13);

	// kryleig return values 
	Eigen::VectorXd d;
	Eigen::MatrixXd V;
	Eigen::VectorXd ev(6);

	for (int i=0; i<N; ++i)
	{
		const int m = i+6;
		std::cout << "subspace dimension m = " << m << std::endl;
		kryleig(A,m, V, d);
		std::sort(d.data(),d.data()+d.size());	   
		ev.head(3) = d.head(3);
		ev.tail(3) = d.tail(3);

		errev(i,0) = m;
		errev.block(i,1,1,6) = ev.transpose();
		errev.block(i,7,1,6) = (ev-ev0).cwiseAbs().transpose();
	}

	{
		mgl::Figure fig;
		fig.plot(errev.col(0), errev.col(1), "r+");
		fig.plot(errev.col(0), errev.col(2), "m+");
		fig.plot(errev.col(0), errev.col(3), "b+");
		fig.title("Krylov subspace iteration - smallest ritz values");
		fig.xlabel("Dimension m of Krylov space");
		fig.ylabel("Ritz value");
		fig.setFontSize(3);
		fig.save("kryleigmin");
	}

	{
		mgl::Figure fig;
		fig.plot(errev.col(0), errev.col(6), "r^");
		fig.plot(errev.col(0), errev.col(5), "m^");
		fig.plot(errev.col(0), errev.col(4), "b^");
		fig.title("Krylov subspace iteration - largest ritz values");
		fig.xlabel("Dimension m of Krylov space");
		fig.ylabel("Ritz value");
		fig.setFontSize(3);
		fig.save("kryleigmax");
	}

	{
		mgl::Figure fig;
		fig.setlog(false, true);
		fig.plot(errev.col(0), errev.col(7), "r+");
		fig.plot(errev.col(0), errev.col(8), "m+");
		fig.plot(errev.col(0), errev.col(9), "b+");
		fig.title("Krylov subspace iteration - error of smallest ritz values");
		fig.xlabel("Dimension m of Krylov space");
		fig.ylabel("Error of Ritz value");
		fig.setFontSize(3);
		fig.save("kryleigminerr");
	}

	{
		mgl::Figure fig;
		fig.setlog(false, true);
		fig.plot(errev.col(0), errev.col(10), "r^");
		fig.plot(errev.col(0), errev.col(11), "m^");
		fig.plot(errev.col(0), errev.col(12), "b^");
		fig.title("Krylov subspace iteration - error of largest ritz values");
		fig.xlabel("Dimension m of Krylov space");
		fig.ylabel("Error of Ritz value");
		fig.setFontSize(3);
		fig.save("kryleigmaxerr");
	}

}
