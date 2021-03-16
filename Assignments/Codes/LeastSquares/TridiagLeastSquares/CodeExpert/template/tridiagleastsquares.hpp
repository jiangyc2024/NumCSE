// TO DO: Implement lsqEst().
// Note: Your code will not compile until after you have defined the function.
// Note: The test anticipates the output in the order (alpha, beta), and not
// (beta, alpha).
#include <Eigen/Dense>
using namespace Eigen;

/* @brief Minimize  |Tz - c|_2 wrt (alpha,beta), where T is triagonal,
 * T.diagonal()=alpha, T.diagonal(+/-1)=beta.
 * @param[in] z An $n$-dimensional vector of data for lhs of Tz=c.
 * @param[in] c An $n$-dimensional vector of data for rhs of Tz=c.
 * @param[out] x The parameters of T, x=(alpha,beta), that solve the least
 * squares problem Tz=c.
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd lsqEst(const VectorXd &z, const VectorXd &c) {
  // START
  return VectorXd();
  // END
}
/* SAM_LISTING_END_0 */
