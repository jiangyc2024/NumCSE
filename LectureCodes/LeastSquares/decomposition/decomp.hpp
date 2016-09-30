# include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
# include <Eigen/QR>

// Computation of \com{full QR-decomposition} \eqref{eq:fullqr},
// dense matrices built for both QR-factors (expensive!)
std::pair<MatrixXd,MatrixXd> qr_decomp_full(const MatrixXd& A) {
  Eigen::HouseholderQR<MatrixXd> qr(A);
  MatrixXd Q = qr.householderQ();
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

# include <Eigen/SVD>

// get U,S,V such that $A = U S V^H$
void svd_decomp(const MatrixXd& A, MatrixXd& U, MatrixXd& S, MatrixXd& V) {
  // Important! If replacing 'ComputeFull' by 'ComputeThin', 
  // the methods matrixU() and matrixV() will return the 
  // matrices U and V but not extended to an orthonormal basis!
  Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = svd.matrixU(); // get complete matrix U 
  V = svd.matrixV(); // get complete matrix V 
  VectorXd sv = svd.singularValues(); // get singular values as vector
  S = MatrixXd::Zero(A.rows(), A.cols());
  const unsigned n = sv.size(); // no. of singular values
  S.block(0,0,n,n) = sv.asDiagonal(); // build matrix S
}
