# include <cmath>
# include <Eigen/Dense>
using Eigen::VectorXd;

// Computation of coefficients \Blue{\alpha, \beta} from \ref{thm:onp:rec}
// IN : \texttt{t} = points in the definition of the discrete $L^2$-inner product
//      \texttt{n} = maximal index desired
// OUT: \texttt{alpha}, \texttŧ{beta} = coefficients of recursion 
void coeffortho(const VectorXd& t, const int n, VectorXd& alpha, VectorXd& beta) {
  const int m = t.size(); // maximal degree of orthogonal polynomial
  alpha = VectorXd( std::min(n-1, m-2) + 1 );
  beta = VectorXd( std::min(n-1, m-2) + 1 );

  alpha(0) = t.sum()/m;

  // initialization of recursion; we store only the values of 
  // the polynomials a the points in \Blue{$\Ct$}
  VectorXd p0,
           p1 = VectorXd::Ones(m),
           p2 = t - alpha(0)*VectorXd::Ones(m);

  // Main loop
  for (unsigned k = 0; k < std::min(n-1, m-2); ++k) {
    p0 = p1; p1 = p2;
    // 3-term recursion \eqref{eq:onp:rec}
    alpha(k + 1) = p1.dot( t.cwiseProduct(p1) ) / std::pow( p1.norm(), 2 );
    beta(k) = std::pow( ( p1.norm() / p0.norm() ), 2 );
    p2 = (t - alpha(k + 1)*VectorXd::Ones(m)).cwiseProduct( p1 ) - beta(k)*p0;
  }
}
