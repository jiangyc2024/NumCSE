#include <complex>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

/* SAM_LISTING_BEGIN_0 */
template <typename Scalar>
void fft2(Eigen::MatrixXcd &C, const Eigen::MatrixBase<Scalar> &Y) {
  using idx_t = Eigen::MatrixXcd::Index;
  const idx_t m = Y.rows(),n=Y.cols();
  C.resize(m,n);
  Eigen::MatrixXcd tmp(m,n);

  Eigen::FFT<double> fft; // Helper class for DFT
  // Transform rows of matrix \Blue{$\VY$}
  for (idx_t k=0;k<m;k++) {
    Eigen::VectorXcd tv(Y.row(k));
    tmp.row(k) = fft.fwd(tv).transpose();
  }

  // Transform columns of temporary matrix
  for (idx_t k=0;k<n;k++) {
    Eigen::VectorXcd tv(tmp.col(k));
    C.col(k) = fft.fwd(tv);
  }
}
/* SAM_LISTING_END_0 */  

/* SAM_LISTING_BEGIN_1 */
template <typename Scalar>
void ifft2(Eigen::MatrixXcd &C, const Eigen::MatrixBase<Scalar> &Y) {
  using idx_t = Eigen::MatrixXcd::Index;
  const idx_t m = Y.rows(),n=Y.cols();
  fft2(C,Y.conjugate()); C = C.conjugate()/(m*n);
}
/* SAM_LISTING_END_1 */  
