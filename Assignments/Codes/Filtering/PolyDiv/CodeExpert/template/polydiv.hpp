#include <iostream>
#include <Eigen/Dense>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <unsupported/Eigen/FFT>

using namespace Eigen;

/* @brief Polynomial multiplication -- naive implementation
 * @param[in] u Vector of coefficients of polynomial $u$
 * @param[in] v Vector of coefficients of polynomial $v$
 * @param[out] uv Vector of coefficients of polynomial $uv = u*v$
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd polyMult_naive(const VectorXd &u, const VectorXd &v) {
  // Fetch degrees of input polynomials
  int degu = u.size() - 1;
  int degv = v.size() - 1;
  // Object for product polynomial p = uv
  int degp = degu + degv;

  VectorXd uv(degp + 1);

  // TODO: multiply polynomials $u$ and $v$ naively
  // START

  // END

  return uv;
}
/* SAM_LISTING_END_0 */

/* @brief Polynomial multiplication -- efficient implementation
 * @param[in] u Vector of coefficients of polynomial $u$
 * @param[in] v Vector of coefficients of polynomial $v$
 * @param[out] uv Vector of coefficients of polynomial $uv = u*v$
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd polyMult_fast(const VectorXd &u, const VectorXd &v) {
  // Initialization
  int m = u.size();
  int n = v.size();

  VectorXd u_tmp = u;
  u_tmp.conservativeResizeLike(VectorXd::Zero(u.size() + n - 1));
  VectorXd v_tmp = v;
  v_tmp.conservativeResizeLike(VectorXd::Zero(v.size() + m - 1));

  VectorXd uv;

  // TODO: multiply polynomials $u$ and $v$ efficiently
  // START

  // END

  return uv;
}
/* SAM_LISTING_END_1 */

/* @brief Polynomial division -- efficient implementation
 * @param[in] uv Vector of coefficients of polynomial $uv$
 * @param[in] u Vector of coefficients of polynomial $u$
 * @param[out] uv Vector of coefficients of polynomial $v$
 * from $uv = u*v$ (division without remainder)
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd polyDiv(const VectorXd &uv, const VectorXd &u) {
  // Initialization
  int mn = uv.size();
  int m = u.size();
  // need well behaved input
  assert(uv(mn - 1) != 0. && u(m - 1) != 0.);
  if (mn < m) {
    std::cerr << "uv can't be divided by u\n";
    return Eigen::VectorXd(0);
  }
  int dim = mn;
  int n = mn - m + 1;

  VectorXd v;

  // TODO: divide polynomials $uv$ and $u$ efficiently
  // START
  
  // END
  return v;
}
/* SAM_LISTING_END_2 */
