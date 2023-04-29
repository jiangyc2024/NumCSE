///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <Eigen/Dense>
# include <cmath>

namespace coeffortho {


using Eigen::VectorXd;
using Eigen::Index;

inline
/* SAM_LISTING_BEGIN_0 */
// Computation of coefficients \Blue{$\alpha$, $\beta$} from \ref{thm:onp:rec}
// IN : \texttt{t} = points in the definition of the discrete $L^2$-inner product
//      \texttt{n} = maximal index desired
//      \texttt{alpha}, \texttt{beta} are used to save coefficients of recursion 
void coeffortho(const VectorXd& t, const Index n, VectorXd& alpha, VectorXd& beta) {
  const Index m = t.size(); // maximal degree of orthogonal polynomial
  alpha = VectorXd( std::min(n-1, m-2) + 1 );
  beta = VectorXd( std::min(n-1, m-2) + 1 );
  alpha(0) = t.sum()/static_cast<double>(m);
  // initialization of recursion; we store only the values of 
  // the polynomials a the points in \Blue{$\Ct$}
  VectorXd p0;
  VectorXd p1 = VectorXd::Ones(m);
  VectorXd p2 = t - alpha(0)*VectorXd::Ones(m);
  for (Index k = 0; k < std::min(n-1, m-2); ++k) {
    p0 = p1; p1 = p2;
    // 3-term recursion \eqref{eq:onp:rec}
    alpha(k+1) = p1.dot(t.cwiseProduct(p1))/p1.squaredNorm();
    beta(k) = p1.squaredNorm()/p0.squaredNorm();
    p2 = (t-alpha(k+1)*VectorXd::Ones(m)).cwiseProduct(p1)-beta(k)*p0;
  }
}
/* SAM_LISTING_END_0 */


} //namespace coeffortho