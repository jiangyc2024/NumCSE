# include <Eigen/Dense>
# include "polyfit.hpp"
# include "polyval.hpp"

using Eigen::VectorXd;
/* Evaluation of the interpolation polynomials with polyfit+polyval
 * IN:  t = nodes
 *      y = values in t
 *      x = evaluation points
 * OUT: v = values of interpolant in x                              */
void ipoleval(const VectorXd& t, const VectorXd& y, const VectorXd& x, VectorXd& v) {
  // get coefficients using polyfit
  VectorXd p =  polyfit(t, y, y.size() - 1);
  // evaluate using polyval (<-> horner scheme)
  polyval(p, x, v);
}
