#include "divdiff.hpp"

/* SAM_LISTING_BEGIN_1 */
// Evaluation of a polynomial in Newton form, that is, represented through the
// vector of its basis expansion coefficients with respect to the Newton basis
// \eqref{eq:newtbas}.
Eigen::VectorXd evalNewtonForm(const Eigen::VectorXd &t,
                               const Eigen::VectorXd &a,
                               const Eigen::VectorXd &x) {
  const unsigned int n = a.size() - 1;
  const Eigen::VectorXd ones = VectorXd::Ones(x.size());
  Eigen::VectorXd p{a[n] * ones};
  for (int j = n - 1; j >= 0; --j) {
    p = (x - t[j] * ones).cwiseProduct(p) + a[j] * ones;
  }
  return p;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_0 */
// Evaluation of polynomial in Newton basis (divided differences)
// IN:  t = nodes (mutually different)
//      y = values in t
//      x = evaluation points (as Eigen::Vector)
// OUT: p = values in x                                           */
void evaldivdiff(const Eigen::VectorXd &t, const Eigen::VectorXd &y,
                 const Eigen::VectorXd &x, Eigen::VectorXd &p) {
  // Get Newton coefficients of polynomial (non in-situ implementation!)
  Eigen::VectorXd coeffs;
  divdiff(t, y, coeffs);
  // evaluate
  p = evalNewtonForm(t, coeffs, x);
}
