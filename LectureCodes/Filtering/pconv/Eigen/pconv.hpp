#include <Eigen/Dense>
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXcd pconv(const Eigen::VectorXcd &u, const Eigen::VectorXcd &x) {
  const int n = x.size();
  Eigen::VectorXcd z = Eigen::VectorXcd::Zero(n);
  // ``naive'' two-loop implementation of discrete periodic convolution
  for (int k = 0; k < n; ++k) {
    for (int j = 0, l = k; j <= k; ++j, --l)  z[k] += u[l] * x[j];
    for (int j = k + 1, l = n - 1; j < n; ++j, --l)  z[k] += u[l] * x[j];
  }
  return z;
}
/* SAM_LISTING_END_0 */
