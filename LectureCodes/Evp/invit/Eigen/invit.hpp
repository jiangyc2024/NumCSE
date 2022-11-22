#include <Eigen/Dense>

namespace invit {


// inverse iteration for computing \Blue{$\lambda_{\min}(\VA)$} and associated eigenvector
inline std::pair<double, Eigen::VectorXd> invit(const Eigen::MatrixXd &A, double tol = 1e-6)
{
	auto lu = A.lu(); // Magenta{single} intial LU-factorization, see Rem.~\ref{rem:seqsolvelse}

	const Eigen::Index n = A.rows();
	Eigen::VectorXd x = Eigen::VectorXd::Random(n); // random initial guess
	x.normalize();

	Eigen::VectorXd y = lu.solve(x);
	double lmin = 1/y.norm();
	y = y*lmin;
	double lold = 0;

	while (std::abs(lmin-lold) > tol*lmin) // termination, if small \emph{relative change}
	{
		lold = lmin; 
		x = y; 
		y = lu.solve(x); // core iteration: \Blue{$\Vy = \VA^{-1}\Vx$}, 
	 	lmin = 1./y.norm(); // new  approxmation of \Blue{$\lambda_{\min}(\VA)$}
		y = y*lmin; // normalization \Blue{$\Vy := \frac{\Vy}{\N{\Vy}_2}$}
	} 

	return std::make_pair(lmin, y);
}


} //namespace invit
