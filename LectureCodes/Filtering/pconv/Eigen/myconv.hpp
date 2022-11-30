#include "pconvfft.hpp"

#include <algorithm>
#include <cassert>

inline
// Discrete convolution of two vectors of equal length by periodic convolution
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXcd fastconv(const Eigen::VectorXcd &h,
                          const Eigen::VectorXcd &x) {
  assert(x.size() == h.size());
  const Eigen::Index n = h.size();
  // Zero padding, cf. \eqref{eq:zeropad}, and periodic discrete convolution
  //  of length \Blue{$2n-1$}, \cref{cpp:pconvfft}
  return pconvfft(
      (Eigen::VectorXcd(2 * n - 1) << h, Eigen::VectorXcd::Zero(n - 1))
          .finished(),
      (Eigen::VectorXcd(2 * n - 1) << x, Eigen::VectorXcd::Zero(n - 1))
          .finished());
}
/* SAM_LISTING_END_0 */

inline
// General discrete convolution: double loop implementation 
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXcd dconv(const Eigen::VectorXcd &h, const Eigen::VectorXcd &x) {
  using idx_t = Eigen::Index;
  const idx_t n = h.size();
  const idx_t m = x.size();
  Eigen::VectorXcd y{Eigen::VectorXcd::Zero(n + m - 1)};
  for (idx_t k = 0; k < n + m - 1; ++k) {
    for (idx_t j = std::max(0L, k - n + 1); j < std::min(m, k + 1); ++j) {
      y[k] += h[k - j] * x[j];
    }
  }
  return y;
}
/* SAM_LISTING_END_1 */
