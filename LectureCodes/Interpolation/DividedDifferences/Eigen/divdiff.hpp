#include <Eigen/Dense>

using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// IN:  t = node set (mutually different)
//      y = nodal values
// OUT: c = coefficients of polynomial in Newton basis
void divdiff(const VectorXd &t, const VectorXd &y, VectorXd &c) {
  const int n = y.size() - 1;
  c = y;
  // Follow scheme \eqref{eq:ddscheme}, recursion \eqref{eq:acrec}
  for (int l = 1; l < n; ++l)
    for (int j = n - l; j < n; ++j)
      c[j] = (c[j] - c[j - 1]) / (t[j] - t[n - 1 - l]);
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
// IN:  t = node set (mutually different)
//      y = nodal values
// OUT: y = coefficients of polynomial in Newton basis
void divdiff(const VectorXd &t, VectorXd &y) {
  const int n = y.size() - 1;
  // Follow scheme \eqref{eq:ddscheme}, recursion \eqref{eq:acrec}
  for (int l = 1; l < n; ++l)
    for (int j = n - l; j < n; ++j)
      y[j] = (y[j] - y[j - 1]) / (t[j] - t[n - 1 - l]);
}
/* SAM_LISTING_END_1 */
