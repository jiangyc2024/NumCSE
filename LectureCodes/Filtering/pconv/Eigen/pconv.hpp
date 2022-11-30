#include <Eigen/Dense>

inline
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXcd pconv(const Eigen::VectorXcd &u, const Eigen::VectorXcd &x) {
  const Eigen::Index n = x.size();
  Eigen::VectorXcd z = Eigen::VectorXcd::Zero(n);
  // ``naive'' two-loop implementation of discrete periodic convolution
  for (Eigen::Index k = 0; k < n; ++k) {
    for (Eigen::Index j = 0, l = k; j <= k; ++j, --l) {
      z[k] += u[l] * x[j];
    }
    for (Eigen::Index j = k + 1, l = n - 1; j < n; ++j, --l) {
      z[k] += u[l] * x[j];
    }
  }
  return z;
}
/* SAM_LISTING_END_0 */
