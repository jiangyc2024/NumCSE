// TO DO: Implement linReg() and expFit().
// Note: The tests anticipate the outputs in the order (alpha, beta), and not (beta, alpha).
#include <Eigen/Dense>
using namespace Eigen;

/* @brief Solve the linear regression problem (fitting a line to data)
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of the line fitted
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd linReg(const VectorXd &t, const VectorXd &y) {
  // START

  // END
}
/* SAM_LISTING_END_0 */

/* @brief Solve the linearized exponential problem
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of the line fitted to the linearized problem
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd expFit(const VectorXd &t, const VectorXd &y) {
  // START

  // END
}
/* SAM_LISTING_END_1 */
