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
  for (int i = 0; i <= degp; ++i) {
    int fst = std::max(0, i - degv);
    int lst = std::min(degu, i);
    for (int j = fst; j <= lst; ++j) {
      uv(i) += u(j) * v(i - j);
    }
  }
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
  Eigen::FFT<double> fft;
  VectorXcd u_tmp_ = fft.fwd(u_tmp);
  VectorXcd v_tmp_ = fft.fwd(v_tmp);
  VectorXcd tmp = u_tmp_.cwiseProduct(v_tmp_);
  uv = fft.inv(tmp).real();
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

  // zero padding
  VectorXd uv_tmp = uv;
  uv_tmp.conservativeResizeLike(VectorXd::Zero(dim));
  VectorXd u_tmp = u;
  u_tmp.conservativeResizeLike(VectorXd::Zero(dim));

  Eigen::FFT<double> fft;
  VectorXcd uv_tmp_ = fft.fwd(uv_tmp);
  VectorXcd u_tmp_ = fft.fwd(u_tmp);

  // check divisibility: case(i)
  for (int i = 0; i < dim; ++i) {
    if (abs(uv_tmp_(i)) < 1e-13) {
      if (abs(u_tmp_(i)) < 1e-13) {
        // make cwiseQuotient at i-th position equal 0
        uv_tmp_(i) = 0.; // complex assignment (0., 0.)
        u_tmp_(i) = 1.;  // complex assignment (1., 0.)
      } else {
        std::cerr << "uv can't be divided by u\n";
        return Eigen::VectorXd(0);
      }
    }
  }

  VectorXcd tmp = uv_tmp_.cwiseQuotient(u_tmp_);

  v = fft.inv(tmp).real();
  // check divisibility: case(ii)
  for (int i = n; i < dim; ++i) {
    if (abs(v(i)) > 1e-13) {
      std::cerr << "uv can't be divided by u\n";
      return Eigen::VectorXd(0);
    }
  }

  // reshape v to a suitable size
  v.conservativeResizeLike(VectorXd::Zero(n));
  // (mn-1) - (m-1) + 1
  // END
  return v;
}
/* SAM_LISTING_END_2 */
