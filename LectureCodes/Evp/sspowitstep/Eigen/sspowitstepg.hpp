#include <Eigen/Dense>

/* power iteration with orthonormalization for \Blue{$\VA=\VA^T$}.
 * columns of matrix \Blue{$\VV$} span subspace for power iteration.
 * V = A*V;         % actual power iteration on individual columns
 * [V,R] = qr(V,0); % Gram-Schmidt orthonormalization \eqref{eq:GSO} */
inline Eigen::MatrixXd sspowitstepg(const Eigen::MatrixXd &A, const Eigen::MatrixXd &V)
{
	const Eigen::MatrixXd B = A*V;
	auto qr = B.householderQr();
	return qr.householderQ();
}
