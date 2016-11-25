#include <Eigen/Dense>

/* Note: this implementation of the Lanczos process also records the orthonormal CG residuals in the columns of the matrix \texttt{V}, which is not needed when only eigenvalue approximations are desired. */
void lanczos(const Eigen::MatrixXd& A, int k, const Eigen::VectorXd &z0, Eigen::MatrixXd &V, Eigen::VectorXd& alpha, Eigen::VectorXd& beta)
{
	V.resize(A.rows(),k+1);
	alpha.resize(k);
	beta.resize(k);

	V.col(0) = z0 / z0.norm();

	// Vectors storing entries of tridiagonal matrix \eqref{eq:tdcg}
	alpha = Eigen::VectorXd::Zero(k,1); 
	beta = Eigen::VectorXd::Zero(k,1);

	for (int j=0; j<k; ++j)
	{
		Eigen::VectorXd q = A*V.col(j);
		alpha(j) = q.dot(V.col(j));
		Eigen::VectorXd w = q - alpha(j)*V.col(j);

		if (j > 0) w = w - beta(j-1)*V.col(j-1);

		beta(j) = w.norm();
		V.col(j+1) = w/beta(j);
	}
	beta = beta.head(k-1);
}
