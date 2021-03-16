#include "fft2.hpp" // contains our implementation of Matlab's fft2
#include <Eigen/Dense>
#include <iostream>
using Eigen::MatrixXcd;
using Eigen::MatrixXd;
// typedef to avoid writing the whole type
typedef std::complex<double> complex;

/* SAM_LISTING_BEGIN_0 */
MatrixXd deblur(const MatrixXd &C, const MatrixXd &S, const double tol = 1e-3) {
  const long m = C.rows(), n = C.cols(), M = S.rows(), N = S.cols();
  const long L = (M - 1) / 2;
  if (M != N) {
    throw std::runtime_error("Error: S not quadratic!");
  }
  MatrixXd Spad = MatrixXd::Zero(m, n);
  // Zero padding, see \eqref{eq:Mdeblur}
  Spad.block(0, 0, L + 1, L + 1) = S.block(L, L, L + 1, L + 1);
  Spad.block(m - L, n - L, L, L) = S.block(0, 0, L, L);
  Spad.block(0, n - L, L + 1, L) = S.block(L, 0, L + 1, L);
  Spad.block(m - L, 0, L, L + 1) = S.block(0, L, L, L + 1);
  // Inverse of blurring operator (fft2 expects a complex matrix)
  MatrixXcd SF = fft2(Spad.cast<complex>());
  // Test for invertibility
  if (SF.cwiseAbs().minCoeff() < tol * SF.cwiseAbs().maxCoeff()) {
    std::cerr << "Error: Deblurring impossible!\n";
  }
  // DFT based deblurring
  return fft2(ifft2(C.cast<complex>()).cwiseQuotient(SF)).real();
}
/* SAM_LISTING_END_0 */
