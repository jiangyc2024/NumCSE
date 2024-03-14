///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>

/* SAM_LISTING_BEGIN_0 */
// Generic numerical quadrature routine implementing \eqref{eq:quadform}: 
// \texttt{f} is a handle to a function, e.g. as lambda function
// \texttt{c}, \texttt{w} pass quadrature nodes \Blue{$\qn_j\in [a,b]$}, and weights \Blue{$\qw_j\in\bbR$}
// in a Eigen::VectorXd
template <class Function>
double quadformula(Function &&f, const Eigen::VectorXd& c,const Eigen::VectorXd& w) {
  const Eigen::Index n = c.size();
  double I = 0;
  for (Eigen::Index i = 0; i < n; ++i) { I += w(i)*f(c(i)); }
  return I;
}
/* SAM_LISTING_END_0 */
