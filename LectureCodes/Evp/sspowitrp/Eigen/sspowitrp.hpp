#include <Eigen/Dense>
#include <eigensolversort.hpp>

/* Power iteration with Ritz projection for matrix \Blue{$\VA=\VA^T\in\bbR^{n,n}$}:
 * Subspace of dimension \Blue{$m\leq n$} is used to compute the \Blue{$k\leq m$} largest
 * eigenvalues of \Blue{$\VA$} and associated eigenvectors.
 */
void sspowitrp(const Eigen::MatrixXd &A, int k, int m, double tol, int maxIt, Eigen::VectorXd &maxEv, Eigen::MatrixXd &maxV)
{
	assert(k <= m);

	int n = A.rows();

	Eigen::MatrixXd V = Eigen::MatrixXd::Identity(n,m); // eigenvectors
	Eigen::VectorXd ev; // eigenvalues

	// (Arbitrary) initial eigenvectors
	Eigen::VectorXd d = Eigen::VectorXd(m);

	// The approximate eigenvectors are stored in the columns of \Blue{$\VV\in\bbR^{n,m}$}
	for (int i=0; i<maxIt; ++i)
	{
		// Power iteration and orthonormalization
		auto qr = (A*V).householderQr();
		Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n,m);

		// Small \Blue{$m\times m$} eigenvalue problem for Ritz projection
		Eigen::EigenSolver<Eigen::MatrixXd> esolver(Q.transpose()*A*Q);

		// recover approximate eigenvectors
		Eigen::MatrixXd U; 
		std::tie(ev, U) = eigensolversort(esolver, false);

		// 2nd part of Ritz projection
		V = Q*U;

		std::cout << "err="<< (ev-d).norm() << std::endl;
		if ((ev-d).norm() < tol*ev.cwiseAbs().maxCoeff()) break;

		d = ev;
	}

	maxEv.resize(k);
	maxV.resize(n,k);

	maxEv = ev.tail(k);
	maxV = V.rightCols(k); 
}

