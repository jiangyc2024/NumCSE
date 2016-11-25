#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <eigensolversort.hpp>

Eigen::MatrixXd lanczosev(const Eigen::MatrixXd& A, int k, int j, const Eigen::VectorXd& w)
{
	Eigen::MatrixXd V(A.rows(), k+1);
	V.col(0) = w / w.norm();

	Eigen::VectorXd alpha = Eigen::VectorXd::Zero(k);
	Eigen::VectorXd beta = Eigen::VectorXd::Zero(k);

	Eigen::MatrixXd result(k, j+1);
	Eigen::VectorXd vt;

	for (int l=0; l<k; ++l)
	{
		vt = A*V.col(l);
		if (l > 0) vt = vt - beta(l-1)*V.col(l-1);

		alpha(l) = vt.dot(V.col(l));
		vt = vt - alpha(l)*V.col(l);
		beta(l) = vt.norm();
		V.col(l+1) = vt/beta(l);

		Eigen::MatrixXd T2 = Eigen::MatrixXd::Zero(l+1,l+1);

		// create diagonals
		if (l>1)
		{
			Eigen::VectorXd dl(l);
			Eigen::VectorXd du(l);

			du = beta.head(l);
			dl = beta.head(l);

			T2.diagonal(-1) = dl;
			T2.diagonal(1) = du;
		}

		T2.diagonal(0) = alpha.head(l+1);

		// eigendecomposition
		Eigen::EigenSolver<Eigen::MatrixXd> esolver(T2);
		Eigen::MatrixXd U; 
		Eigen::VectorXd ev; 
		std::tie(ev, U) = eigensolversort(esolver);

		Eigen::VectorXd ev2(j);

		if (j > l)
		{
			ev2.head(j-l-1) = Eigen::VectorXd::Zero(j-l-1);
			ev2.tail(l+1) = ev;
		}
		else
		{
			ev2 = ev.tail(j);
		}
		result(l,0) = l;
		result.row(l).tail(j) = ev2;
	}

	return result;
}
