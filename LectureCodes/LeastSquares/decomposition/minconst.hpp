///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <limits>
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;

# include <Eigen/SVD>
const double EPS = std::numeric_limits<double>::epsilon();

/* SAM_LISTING_BEGIN_0 */
// \eigen based function for solving \eqref{lsq:minconst};
// minimizer returned nin x, mininum as return value
double minconst(VectorXd &x,const MatrixXd &A) {
  MatrixXd::Index m=A.rows(),n=A.cols();
  if (m < n) throw std::runtime_error("A must be tall matrix");
  // SVD factor \Blue{$\VU$} is \textbf{not} computed!
  Eigen::JacobiSVD<MatrixXd> svd(A,Eigen::ComputeThinV);
  x.resize(n); x.setZero(); x(n-1) = 1.0; // \Blue{$\Ve_n$}
  x = svd.matrixV()*x;
  return (svd.singularValues())(n-1); 
}
/* SAM_LISTING_END_0 */

