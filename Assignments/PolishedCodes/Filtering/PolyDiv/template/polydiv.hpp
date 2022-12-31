#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unsupported/Eigen/FFT>

/**
 * @brief Polynomial multiplication -- naive implementation
 *
 * @param u Vector of coefficients of polynomial $u$
 * @param v Vector of coefficients of polynomial $v$
 * @return uv Vector of coefficients of polynomial $uv = u*v$
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd polyMult_naive(const Eigen::VectorXd &u,
                               const Eigen::VectorXd &v) {
  // Fetch degrees of input polynomials
  const int degu = u.size() - 1;
  const int degv = v.size() - 1;

  // Object for product polynomial p = uv
  const int degp = degu + degv;

  Eigen::VectorXd uv(degp + 1);

  // TODO: (4-3.a) Multiply polynomials $u$ and $v$ naively.
  // START

  // END

  return uv;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Polynomial multiplication -- efficient implementation
 *
 * @param u Vector of coefficients of polynomial $u$
 * @param v Vector of coefficients of polynomial $v$
 * @return uv Vector of coefficients of polynomial $uv = u*v$
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd polyMult_fast(const Eigen::VectorXd &u,
                              const Eigen::VectorXd &v) {
  // Initialization
  const unsigned int m = u.size();
  const unsigned int n = v.size();

  Eigen::VectorXd u_tmp(m + n - 1);
  u_tmp.head(m) = u;
  u_tmp.tail(n - 1).setZero();

  Eigen::VectorXd v_tmp(m + n - 1);
  v_tmp.head(n) = v;
  v_tmp.tail(m - 1).setZero();

  Eigen::VectorXd uv;

  // TODO: (4-3.b) Multiply polynomials $u$ and $v$ efficiently.
  // START

  // END

  return uv;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Polynomial division -- efficient implementation
 *
 * @param uv Vector of coefficients of polynomial $uv$
 * @param u Vector of coefficients of polynomial $u$
 * @return uv Vector of coefficients of polynomial $v$
 * from $uv = u*v$ (division without remainder)
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd polyDiv(const Eigen::VectorXd &uv, const Eigen::VectorXd &u) {
  // Initialization
  const unsigned int mn = uv.size();
  const unsigned int m = u.size();
  // need well behaved input
  if (mn < m) {
    std::cerr << "uv can't be divided by u\n";
    return Eigen::VectorXd(0);
  }
  const unsigned int dim = mn;
  const unsigned int n = mn - m + 1;

  Eigen::VectorXd v;

  // TODO: (4-3.e) Divide polynomials $uv$ and $u$ efficiently.
  // START

  // END

  return v;
}
/* SAM_LISTING_END_2 */
