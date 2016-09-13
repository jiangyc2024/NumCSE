# include <Eigen/Dense>
# include <Eigen/SVD>
# include <Eigen/QR>

using Eigen::VectorXd;
using Eigen::MatrixXd;

// get Q,R such that $A = Q R$
void qr_decomp(const MatrixXd& A, MatrixXd& Q, MatrixXd& R) {
  Eigen::HouseholderQR<MatrixXd> qr(A);
  Q = qr.householderQ();
  R = qr.matrixQR().template triangularView<Eigen::Upper>();
}

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
  const unsigned n = s.size(); // no. of singular values
  S.block(0,0,n,n) = sv.asDiagonal(); // build matrix S
}
