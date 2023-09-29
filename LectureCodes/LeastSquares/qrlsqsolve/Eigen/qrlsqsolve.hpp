///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <Eigen/QR>


namespace qrlsqsolve {


using Eigen::VectorXd;
using Eigen::MatrixXd;

inline
/* SAM_LISTING_BEGIN_0 */
// Solution of linear least squares problem \eqref{eq:LSQ1} by means of QR-decomposition
// Note: \Blue{$\VA\in\bbR^{m,n}$} with \Blue{$m>n$}, \Blue{$\operatorname{rank}(\VA) = n$} is assumed
// Least squares solution returned in \texttt{x}, residual norm as return value
double qrlsqsolve(const MatrixXd& A, const VectorXd& b,
		  VectorXd& x) {
  const unsigned m = A.rows();
  const unsigned n = A.cols();

  MatrixXd Ab(m, n + 1); Ab << A, b; // Form extended matrix \Blue{$[\VA,\Vb]$} \Label[line]{qrl:0}
  
  // QR-decomposition of extended matrix automatically transforms \Blue{$\Vb$}
  MatrixXd R = Ab.householderQr().matrixQR().template
    triangularView<Eigen::Upper>(); // \Label[line]{qrl:1}

  MatrixXd R_nn = R.block(0, 0, n, n); // R-factor \Blue{$\VR_0$}
  // Compute least squares \com{solution} \Blue{${\Vx} = (\VR)_{1:n,1:n}^{-1}(\VQ^{\top}\Vb)_{1:n}$}
  x = R_nn.template triangularView<Eigen::Upper>().solve(R.block(0, n, n, 1));     
  return R(n, n); // \com{residual norm} \Blue{$ = \N{\VA\wh{\Vx}-\Vb}_2$} (why ?) \label{qrl:2}
}
/* SAM_LISTING_END_0 */

inline
/* SAM_LISTING_BEGIN_1 */
// Solving a full-rank least squares problem \Blue{$\N{\VA\Vx-\Vb}_2\to\min$} in \eigen
double lsqsolve_eigen(const MatrixXd& A, const VectorXd& b,
		      VectorXd& x) {
  x = A.householderQr().solve(b);
  return ((A*x-b).norm());
}
/* SAM_LISTING_END_1 */


// for the keyword "template", see http://eigen.tuxfamily.org/dox/TopicTemplateKeyword.html


} // namespace qrlsqsolve