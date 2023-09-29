#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>

namespace lsqtotal {


using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::cerr;

inline
/* SAM_LISTING_BEGIN_0 */
// computes only solution \Blue{$\Vx$} of fitted consistent LSE
VectorXd lsqtotal(const MatrixXd& A, const VectorXd& b) {
  const unsigned m = A.rows();
  const unsigned n = A.cols(); 
  MatrixXd C(m, n + 1);  C << A, b; // \Blue{$\VC = [\VA,\Vb]$}
  // We need only the SVD-factor \Blue{$\VV$}, see \eqref{tlsq:1}
  MatrixXd V = C.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).matrixV(); 

  // Compute solution according to \eqref{tlsq:3};
  const double s = V(n, n);
  if (std::abs(s) < 1.0E-15) { cerr << "No solution!\n"; return {}; }
  return (-V.col(n).head(n) / s); 
}
/* SAM_LISTING_END_0 */


} // namespace lsqtotal