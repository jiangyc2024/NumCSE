///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// returns the values of the first n - 1 legendre polynomials
// in point x as columns of the matrix L 
void legendre(const unsigned n, const VectorXd& x, MatrixXd& L) {
  L = MatrixXd::Ones(n,n); // \Blue{$p_0(x) = 1$}
  L.col(1) = x;            // \Blue{$p_1(x) = x$}
  for (unsigned j = 1; j < n - 1; ++j) {
    // \Blue{$p_{j+1}(x) = \frac{2j + 1}{j + 1} x p_{j}(x) - \frac{j}{j + 1} p_{j - 1}(x)$} \cref{eq:Legpol}
    L.col(j+1) = (2.*j+1)/(j+1.)*L.col(j-1).cwiseProduct(x)
      -j/(j+1.)*L.col(j-1);
  }
}
/* SAM_LISTING_END_0 */
