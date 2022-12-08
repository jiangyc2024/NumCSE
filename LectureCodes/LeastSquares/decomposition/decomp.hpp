///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
# include <Eigen/QR>

// Computation of \com{full QR-decomposition} \eqref{eq:fullqr},
// dense matrices built for both QR-factors (expensive!)
std::pair<MatrixXd,MatrixXd> qr_decomp_full(const MatrixXd& A) {
  Eigen::HouseholderQR<MatrixXd> qr(A);
  MatrixXd Q = qr.householderQ(); // \Label[line]{cqr:f1}
  MatrixXd R = qr.matrixQR().template triangularView<Eigen::Upper>();
  return std::pair<MatrixXd,MatrixXd>(Q,R);
}

// Computation of \com{economical QR-decomposition} \eqref{eq:ecoqr},
// dense matrix built for Q-factor (possibly expensive!)
std::pair<MatrixXd,MatrixXd> qr_decomp_eco(const MatrixXd& A) {
  using index_t = MatrixXd::Index;
  const index_t m = A.rows(),n = A.cols();
  Eigen::HouseholderQR<MatrixXd> qr(A); 
  MatrixXd Q = (qr.householderQ()*MatrixXd::Identity(m,n)); // \Label[line]{cqr:e1}
  MatrixXd R = qr.matrixQR().block(0,0,n,n).template triangularView<Eigen::Upper>(); // \Label[line]{cqr:e2}
  return std::pair<MatrixXd,MatrixXd>(Q,R);
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
#include <Eigen/SVD>
using MatrixXd = Eigen::MatrixXd;
// Computation of (full) SVD \Blue{$\VA=\VU\Sigmabf\VV^\herm$} $\to$ \cref{thm:svd}
// SVD factors are returned as dense matrices in natural order
std::tuple<MatrixXd,MatrixXd,MatrixXd> svd_full(const MatrixXd& A) {
  Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
  MatrixXd U = svd.matrixU(); // get unitary (square) matrix \Blue{$\VU$} 
  MatrixXd V = svd.matrixV(); // get unitary (square) matrix \Blue{$\VV$} 
  VectorXd sv = svd.singularValues(); // get singular values as vector
  MatrixXd Sigma = MatrixXd::Zero(A.rows(), A.cols());
  const unsigned p = sv.size(); // no. of singular values
  Sigma.block(0,0,p,p) = sv.asDiagonal(); // set diagonal block of \Blue{$\Sigmabf$}
  return std::tuple<MatrixXd,MatrixXd,MatrixXd>(U,Sigma,V);
}

// Computation of economical (thin) SVD \Blue{$\VA=\VU\Sigmabf\VV^\herm$}, see \eqref{eq:svdeco}
// SVD factors are returned as dense matrices in natural order
std::tuple<MatrixXd,MatrixXd,MatrixXd> svd_eco(const MatrixXd& A) {
  Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  MatrixXd U = svd.matrixU(); // get matrix \Blue{$\VU$} with orthonormal columns
  MatrixXd V = svd.matrixV(); // get  matrix \Blue{$\VV$}  with orthonormal columns
  VectorXd sv = svd.singularValues(); // get singular values as vector
  MatrixXd Sigma = sv.asDiagonal(); // build diagonal matrix \Blue{$\Sigmabf$}
  return std::tuple<MatrixXd,MatrixXd,MatrixXd>(U,Sigma,V);
}
/* SAM_LISTING_END_1 */

// Important! If replacing 'ComputeFull' by 'ComputeThin', 
// the methods matrixU() and matrixV() will return the 
// matrices U and V but not extended to an orthonormal basis!
