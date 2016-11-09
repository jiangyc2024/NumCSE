# include <Eigen/Dense>

using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
void divdiff(const VectorXd& t, const VectorXd& y, VectorXd& c) {
  // IN:  t = node set (mutually different)
  //      y = nodal values 
  // OUT: c = coefficients of polynomial in Newton basis
  
  c = y;
  const unsigned n = y.size() - 1;
  for (unsigned l = 0; l < n; ++l) 
    for (unsigned j = l; j < n; ++j) 
      c(j+1) = ( c(j+1) - c(l) )/( t(j+1) - t(l) );
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
// IN:  t = node set (mutually different)
//      y = nodal values 
// OUT: y = coefficients of polynomial in Newton basis
void divdiff(const VectorXd& t, VectorXd& y) {
  const unsigned n = y.size() - 1;
  // Follow scheme \eqref{eq:ddscheme}, recursion \eqref{eq:acrec}
  for (unsigned l = 0; l < n; ++l) 
    for (unsigned j = l; j < n; ++j) 
      y(j+1) = (y(j+1)-y(l))/(t(j+1)-t(l));
}
/* SAM_LISTING_END_1 */
