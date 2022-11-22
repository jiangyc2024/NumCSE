#include <Eigen/Dense>
#include <iostream>

#include "fft2.hpp"  // contains our implementation of Matlab's fft2
// typedef to avoid writing the whole type
using complex = std::complex<double>;

inline
/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd deblur(const Eigen::MatrixXd &C, const Eigen::MatrixXd &S,
                       const double tol = 1e-3) {
  typedef Eigen::Index index_t;
  auto dimensions = std::make_tuple(C.rows(), C.cols(), S.rows(), S.cols());
  const auto [ m, n, M, N ] = dimensions; 
  const index_t L = (M - 1) / 2;
  if (M != N) {
    throw std::runtime_error("Error: S not quadratic!");
  }
  Eigen::MatrixXd Spad = Eigen::MatrixXd::Zero(m, n);
  // Zero padding, see \eqref{eq:Mdeblur}
  Spad.block(0, 0, L + 1, L + 1) = S.block(L, L, L + 1, L + 1);
  Spad.block(m - L, n - L, L, L) = S.block(0, 0, L, L);
  Spad.block(0, n - L, L + 1, L) = S.block(L, 0, L + 1, L);
  Spad.block(m - L, 0, L, L + 1) = S.block(0, L, L, L + 1);
  // Inverse of blurring operator (fft2 expects a complex matrix)
  const Eigen::MatrixXcd SF = fft2(Spad.cast<complex>());
  // Test for invertibility
  if (SF.cwiseAbs().minCoeff() < tol * SF.cwiseAbs().maxCoeff()) {
    std::cerr << "Error: Deblurring impossible!\n";
  }
  // DFT based deblurring
  return ifft2(fft2(C.cast<complex>()).cwiseQuotient(SF)).real();
}
/* SAM_LISTING_END_0 */
