# include <iostream>
# include <Eigen/Dense>
# include "fft2.hpp" // contains our implementation of Matlab's fft2
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
// typedef to avoid writing the whole type
typedef std::complex<double> complex;

/* SAM_LISTING_BEGIN_0 */
MatrixXcd deblur(const MatrixXd& C, const MatrixXd& S, const double tol=1e-3) {
  const long m = C.rows(), n = C.cols(),
             M = S.rows(), N = S.cols(), L = (M-1)/2;
  if (M != N) {
    std::cerr << "Error: S not quadratic!\n";
  }

  MatrixXd Spad = MatrixXd::Zero(m,n);
  // Zero padding
  Spad.block(0, 0, L+1, L+1) = S.block(L, L, L+1, L+1);
  Spad.block(m-L, n-L, L, L) = S.block(0, 0, L, L);
  Spad.block(0, n-L, L+1, L) = S.block(L, 0, L+1, L);
  Spad.block(m-L, 0, L, L+1) = S.block(0, L, L, L+1);
  // Inverse of blurring operator (fft2 expects a complex matrix)
  MatrixXcd SF = fft2(Spad.cast<complex>());
  // Test for invertibility
  if (SF.cwiseAbs().minCoeff() < tol*SF.cwiseAbs().maxCoeff()) {
    std::cerr << "Error: Deblurring impossible!\n";
  }
  // DFT based deblurring
  MatrixXcd D = fft2(ifft(C.cast<complex>()).cwiseQuotient(SF));
  return std::move(D);
}
/* SAM_LISTING_END_0 */
