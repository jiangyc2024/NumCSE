#include "speye.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace rqui {


// Rayleigh quotient iteration (for \emph{normal} \Blue{$\VA\in\bbR^{n,n}$})
inline std::pair<Eigen::VectorXd, double> rqui(Eigen::SparseMatrix<double> &A, int maxit = 200, double tol = 1e-3)
{
	double alpha = 0; 
	double lmin = NAN;
	const Eigen::Index n = A.rows();

	Eigen::VectorXd z = Eigen::VectorXd::Random(n);
	z.normalize(); // \Blue{$\Vz^{(0)}$}
	const Eigen::SparseMatrix<double> E = speye<double>(n);

	for (int i=0; i<maxit; ++i)
	{
		const Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A-E*alpha);
		z = solver.solve(z); // \Blue{$\Vz^{(k+1)} = (\VA-\rho_{\VA}(\Vz^{(k)})\VI)^{-1}\Vx^{(k)}$}\label{rqui:1}
		z.normalize();
		
		lmin = (A*z).dot(z); // Computation of \Blue{$\rho_{\VA}(\Vz^{(k+1)})$}
		if (std::abs(alpha-lmin) < tol*lmin) /* Desired relative accuracy reached ?*/ { 
			break; 
		}

		alpha = lmin;
	}
	return std::make_pair(z, lmin);
}


} //namespace rqui