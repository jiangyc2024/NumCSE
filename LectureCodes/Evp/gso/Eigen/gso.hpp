#include <Eigen/Dense>

/* \Hyperlink{GSO}{Gram-Schmidt orthonormalization} of the columns of \Blue{$\VV\in\bbR^{n,m}$}, see
 * \eqref{eq:GSO}. The vectors \Blue{$\Vq_1,\ldots,\Vq_m$} are returned as the columns of
 * the \Magenta{\emph{orthogonal}} matrix \Blue{$\VQ$}. */
Eigen::MatrixXd gso(const Eigen::MatrixXd &V)
{
	int m = V.cols();
	Eigen::MatrixXd Q(V.rows(),m);

	// first vector is will only be normalized
	Eigen::VectorXd q = V.col(0);
	q.normalize();
	Q.col(0) = q;

	for (int l=1; l<m; ++l)
	{
		q = V.col(l);

		// orthogonalization
		for (int k=0; k<l; ++k)
		{
			q = q - Q.col(k).dot(V.col(l)) * Q.col(k);
		}

		q.normalize();
		Q.col(l) = q;
	}

	return Q;
}
