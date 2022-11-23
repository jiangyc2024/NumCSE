#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

Eigen::VectorXcd eigqr(const Eigen::MatrixXd& M, double tol=1e-6)
{
	int n = M.rows();
	Eigen::MatrixXcd A = M.cast<std::complex<double>>();

	while (Eigen::MatrixXcd(A.triangularView<Eigen::StrictlyLower>()).norm() > tol*A.norm())
	{
		// shift by eigenvalue of lower right 2$\times$2 block closest to $(\VA)_{n,n}$    
		Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eig(A.bottomRightCorner(2,2), false);
		auto sc = eig.eigenvalues();

		int si; // index of min shifted eigenvalues
		(sc - Eigen::VectorXcd::Constant(sc.rows(), A(n-1,n-1))).cwiseAbs().minCoeff(&si);
		auto shift = sc(si);

		// qr decomposition
		Eigen::FullPivHouseholderQR<Eigen::MatrixXcd> qr(A - shift*Eigen::MatrixXcd::Identity(n,n));
		auto Q = qr.matrixQ().eval();

		A = Q.transpose() * A * Q; 
	}

	return A.diagonal();
}
