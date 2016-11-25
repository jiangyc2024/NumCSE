#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


// Ritz projection onto Krylov subspace. An orthonormal basis of \Blue{$\Kryl{m}{\VA}{\mathbf{1}}$} is assembled into the columns of \Blue{$\VV$}. 
void kryleig(const Eigen::MatrixXd &A, int m, Eigen::MatrixXd &V, Eigen::VectorXd &d)
{
	int n = A.rows();
	V.resize(n,m);
	d.resize(m);
   	V.col(0) = Eigen::VectorXd::LinSpaced(n,1,n);
	V.col(0).normalize();

	for (int l=1; l<m; ++l)
	{
		V.col(l) = A*V.col(l-1);

		auto qr = V.householderQr();
		Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n,l+1); // economy size qr decomposition

		Eigen::EigenSolver<Eigen::MatrixXd> esolver(Q.transpose()*A*Q);

		V.leftCols(l+1) = Q*esolver.eigenvectors().real();
		d.head(l+1) = esolver.eigenvalues().real();
	}
}
