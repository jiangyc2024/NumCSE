#include <Eigen/Dense>
#include <vector>
#include <tuple>



// preconditioned inverse iteration
template <class AFunc, class InvBFunc>
std::tuple<Eigen::VectorXd, double, std::vector<double>> pinvit(AFunc& evalA, int n, InvBFunc& invB, double tol = 1e-12, int maxit=50)
{
	// \texttt{invB} $\hat{=}$ handle to function implementing preconditioner \Blue{$\VB^{-1}$}
	double rho = 0, rhon; 
	std::vector<double> res;

	Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(n, 1, n); // initial guess
	z.normalize(); 

	Eigen::VectorXd v, r;

	for (int i=0; i<maxit; ++i)
	{
		v = evalA(z);
		rhon = v.dot(z); // Rayleigh quotient
		r = v- rhon*z; // residual
		z = z - invB(r); // iteration according to \eqref{eq:pcinvit} 
		z.normalize();
		res.push_back(rhon); // tracking iteration
		
		if (std::abs(rho-rhon) < tol*std::abs(rhon)) break;
		else rho = rhon; 
	}
	double lmin = evalA(z).dot(z);
	res.push_back(lmin);
	
	return std::make_tuple(z, lmin, res);
}
