#include <Eigen/Dense>

Eigen::MatrixXd sspowitsteprp(const Eigen::MatrixXd& A, const Eigen::MatrixXd& V)
{
	// power iteration applied to columns of \Blue{$\VV$}
	Eigen::MatrixXd B = A*V;

	// orthonormalization, see \cref{sec:evporth}
	auto qr = B.householderQr();
	Eigen::MatrixXd Q = qr.householderQ();

	// Solve Ritz projected \Blue{$m\times m$} eigenvalue problem
	Eigen::EigenSolver<Eigen::MatrixXd> ev(Q.transpose()*A*Q);

	// recover approximate eigenvectors
	Eigen::MatrixXd M = ev.eigenvectors().real();
	return Q*M;
}
