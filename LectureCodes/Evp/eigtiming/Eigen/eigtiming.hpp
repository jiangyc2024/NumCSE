// #define EIGEN_USE_MKL_ALL

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <iostream>
#include <figure/figure.hpp>
#include "../../../../Utils/timer.h"

// numerical experiment to determine to complexity of eigenvalue decompositions 
// for different types of matrices
void eigtiming()
{
	size_t N = 50;

	// prepare matrices
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(N,N);
	Eigen::MatrixXd B = A.transpose() * A; // symmetric matrix
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N,N); // tridiagonal matrix
	for (int i=-1;i<=1;++i) C.diagonal(i) = A.diagonal(i);

	std::vector<Eigen::MatrixXd> matrices = {A, B, C};

	// general benchmark function
	auto benchEigenDec = []( const Eigen::MatrixXd M, bool computeEigenvectors)
	{
		Timer timer; timer.start();
		for (int k=0; k<3; ++k)
		{
			Eigen::EigenSolver<Eigen::MatrixXd> solver(M, computeEigenvectors);
			timer.lap();
		}
		return timer.min();
	};

	// perform benchmarks
	std::vector<Eigen::VectorXd> times(matrices.size()*2, Eigen::VectorXd(N/5));
	Eigen::VectorXd sizes(N/5);
	std::cout << "b "<< std::endl;
	for (size_t i=0; i<N/5; i++)
	{
		size_t n = 5+i*5;
		std::cout << n << std::endl;
		sizes(i) = n;
		for (size_t j=0; j<matrices.size(); ++j)
		{
			Eigen::MatrixXd Mn = matrices[j].topLeftCorner(n,n);
			times[2*j](i) = benchEigenDec(Mn, false);
			times[2*j+1](i) = benchEigenDec(Mn, true);
		}
	}


	// plot benchmarks

	mgl::Figure fig;
	fig.plot(sizes, times[0]).label("eigenvalues, random matrix");
	fig.plot(sizes, times[1]).label("eigenvalues & eigenvectors, random matrix");

	fig.plot(sizes, times[2]).label("eigenvalues, symmetric matrix");
	fig.plot(sizes, times[3]).label("eigenvalues & eigenvectors, symmetric matrix");

	fig.plot(sizes, times[4]).label("eigenvalues, tridiagonal matrix");
	fig.plot(sizes, times[5]).label("eigenvalues & eigenvectors, tridiagonal matrix");

	fig.legend(1, 1);
  	fig.setlog(true, true, false);
	fig.title("eig runtimes");
	fig.xlabel("matrix size n");
	fig.ylabel("time [s]");
	fig.setFontSize(5);
	fig.save("eigtimeingall");
	

	// Result: Compared to the matlab eig() function Eigens EigenSolver class does not
	// seem to optimize for special types of matrices
}
