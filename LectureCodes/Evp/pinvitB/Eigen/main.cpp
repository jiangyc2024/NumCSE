#include "pinvitB.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <tuple>
#include <string.h>
#include <vector>
#include <figure/figure.hpp>
#include "spdiags.hpp"

// calculate iterations or errors for pinvit with matrix of size n and given tol and maxit
void benchmark(int n, double tol, int maxit, bool calcError, int &iterations, int& iterationspc, Eigen::VectorXd& err, Eigen::VectorXd& errpc)
{
	// build matrices

	VectorXi diag_no(5); // indices of diagonals
	diag_no << -n/2, -1, 0, 1, n/2;

	Eigen::MatrixXd M(n,5); // matrix of diagonals
	M.col(0) = Eigen::VectorXd::Constant(n, 1./n);
	M.col(1) = Eigen::VectorXd::Constant(n, -1);
	M.col(2) = Eigen::VectorXd::Constant(n, 2*(1+1./n));
	M.col(3) = Eigen::VectorXd::Constant(n, -1);
	M.col(4) = Eigen::VectorXd::Constant(n, 1./n);

	Eigen::SparseMatrix<double> A = spdiags(M, diag_no,n,n);
	Eigen::SparseMatrix<double> B = spdiags(M.middleCols(1,3).eval(), diag_no.middleRows(1,3).eval(), n, n);


	// lambda evaluation for A
	auto evalA = [&](const Eigen::VectorXd& x){return A*x; };

	Eigen::SparseLU<Eigen::SparseMatrix<double>> luA(A);
	auto evalAinv = [&](const Eigen::VectorXd& b){
		return luA.solve(b); 
	};

	// lambda evaluation for B
	auto evalB = [&](const Eigen::VectorXd& x){return B*x;};

	Eigen::SparseLU<Eigen::SparseMatrix<double>> luB(B);
	auto evalBinv = [&](const Eigen::VectorXd& b){
		return luB.solve(b); 
	};


	// return values of pinvit
	double lmin;
	Eigen::VectorXd y;
	std::vector<double> res;

	// benchmark matrix A
	std::tie(y, lmin, res) = pinvit(evalA, n, evalAinv, tol, maxit);
	iterations = res.size();

	if (calcError)
	{
		// compute reference eigenvalue with eigensolver
		Eigen::EigenSolver<Eigen::MatrixXd> solver((Eigen::MatrixXd(A)));
		double lminEx = solver.eigenvalues().real().minCoeff();

		Eigen::VectorXd resVec = Eigen::VectorXd::Map(res.data(), res.size());
		err = (resVec - Eigen::VectorXd::Constant(resVec.rows(), lminEx))/ lminEx;
	}


	// benchmark matrix B
	std::tie(y, lmin, res) = pinvit(evalB, n, evalBinv, tol, maxit);
	iterationspc = res.size();

	if (calcError)
	{
		// compute reference eigenvalue with eigensolver
		Eigen::EigenSolver<Eigen::MatrixXd> solverB((Eigen::MatrixXd(B)));
		double lminEx = solverB.eigenvalues().real().minCoeff();

		Eigen::VectorXd resVec = Eigen::VectorXd::Map(res.data(), res.size());
		errpc = (resVec - Eigen::VectorXd::Constant(resVec.rows(), lminEx))/ lminEx;
	}

}

void benchmarkError()
{
	mgl::Figure fig;
	fig.setlog(false, true);

	int N[] = {50,100,200};
	
	for (int n : N)
	{
		std::cout << "Calculate error for n=" << n << std::endl;
		Eigen::VectorXd err;
		Eigen::VectorXd errpc;

		int nit;
		int nitpc;
		benchmark(n, 1e-14, 50, true, nit, nitpc, err, errpc);

		Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(nit, 0, nit-1);

		std::string label = "INVIT, n = " + std::to_string(n);
		fig.plot(x, err).label(label);

		x = Eigen::VectorXd::LinSpaced(nitpc, 0, nitpc-1);
		std::string label2 = "PINVIT, n = " + std::to_string(n);
		fig.plot(x, errpc).label(label2);
	}

	fig.legend(0, 0.1);
	fig.title("error norms of pinvit");
	fig.xlabel("iteration step");
	fig.ylabel("error in approximation of lambda");
  	fig.setFontSize(5);
	fig.save("pinvitB");
}

void benchmarkIterations()
{
	// number of iterations
	std::vector<int> it;
	std::vector<int> itpc;
	std::vector<int> N; // matrix size
	Eigen::VectorXd err, errpc;

	for (int k=4; k<16; ++k)
	{
		int n = std::pow(2,k);
		N.push_back(n);
		std::cout << "Calculate iterations for n=" << n << std::endl;

		Eigen::VectorXd err;
		Eigen::VectorXd errpc;
		
		int nit, nitpc;
		benchmark(n, 1e-4, 50, false, nit, nitpc, err, errpc);

		// keep track of number of iterations
		it.push_back(nit); 
		itpc.push_back(nitpc);
	}
	mgl::Figure fig;
	fig.setlog(false, true);

	fig.legend(0, 0.1);
	fig.title("number of iterations of pinvit");

	fig.plot(N, it, "b+").label("INVIT");
	fig.plot(N, itpc, "r+").label("PINVIT");

	fig.xlabel("n");
	fig.ylabel("iteration steps");
  	fig.setFontSize(5);
	fig.save("pinvitBitnum");
}

int main()
{
	benchmarkError();
	benchmarkIterations();
}
