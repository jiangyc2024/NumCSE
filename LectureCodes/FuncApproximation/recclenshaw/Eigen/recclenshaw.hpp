///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <Eigen/Dense>
using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// Recursive evaluation of a polynomial \Blue{$p = \sum_{j=1}^{n+1}a_j T_{j-1}$} at point \texttt{x}
// based on \eqref{eq:cstr}
// IN : Vector of coefficients a
//      evaluation point x
// OUT: Value at point x
double recclenshaw(const VectorXd& a, const double x) {
  const VectorXd::Index n = a.size() - 1;
  if      (n == 0) return a(0);            // Constant polynomial
  else if (n == 1) return (x*a(1) + a(0)); // Value \Blue{$\alpha_1*x + \alpha_0$}
  else {
    VectorXd new_a(n);
    new_a << a.head(n - 2), a(n - 2) - a(n), a(n - 1)+ 2*x*a(n);
    return recclenshaw(new_a,x); // recursion 
  }
}
/* SAM_LISTING_END_0 */
