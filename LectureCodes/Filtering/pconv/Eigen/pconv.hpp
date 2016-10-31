# include <Eigen/Dense>
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
VectorXcd pconv(const VectorXcd& u, const VectorXcd& x) {
  using idx_t = VectorXcd::Index; // may be unsigned !
  const idx_t n = x.size(); 
  VectorXcd z = VectorXcd::Zero(n);
  // Need signed indices when differences are formed
  for (long k = 0; k < n; ++k) {
    for (long j = 0; j < n; ++j) {
      long ind = (k - j < 0 ? n + k - j : k - j);
      z(k) += u(ind)*x(j);
    }}
  return z;
}
/* SAM_LISTING_END_0 */
