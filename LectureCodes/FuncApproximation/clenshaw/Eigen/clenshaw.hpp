///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <iostream>
# include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
// Clenshaw algorithm for evaluating \Blue{$p = \sum_{j=1}^{n+1}a_j T_{j-1}$}
// at points passed in vector \texttt{x}
// IN : \Blue{$\mathbf{a}=\left[\alpha_j\right]$}, coefficients for \Blue{$p = \sum_{j=1}^{n+1} \alpha_j T_{j-1}$}
//      \texttt{x} = (many) evaluation points 
// OUT: values $p(x_j)$ for all \Blue{$j$}
VectorXd clenshaw(const VectorXd& a, const VectorXd& x) {
  const int n = a.size() - 1;  // degree of polynomial
  MatrixXd d(n + 1, x.size()); // temporary storage for intermediate values
  for (int c = 0; c < x.size(); ++c) d.col(c) = a; 
  for (int j = n - 1; j > 0; --j) {
    d.row(j) += 2*x.transpose().cwiseProduct(d.row(j+1)); // see \eqref{eq:cstr}
    d.row(j-1) -= d.row(j + 1);
  }
  return d.row(0) + x.transpose().cwiseProduct(d.row(1));
}
/* SAM_LISTING_END_0 */
