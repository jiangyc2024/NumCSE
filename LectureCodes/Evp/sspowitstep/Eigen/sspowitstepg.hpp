#include <Eigen/Dense>

/* power iteration with orthonormalization for \Blue{$\VA=\VA^T$}.
 * columns of matrix \Blue{$\VV$} span subspace for power iteration.
 * V = A*V;         % actual power iteration on individual columns
 * [V,R] = qr(V,0); % Gram-Schmidt orthonormalization \eqref{eq:GSO} */
Eigen::MatrixXd sspowitstepg(const Eigen::MatrixXd &A, const Eigen::MatrixXd &V)
{
	Eigen::MatrixXd B = A*V;
	auto qr = B.householderQr();
	return qr.householderQ();
}
