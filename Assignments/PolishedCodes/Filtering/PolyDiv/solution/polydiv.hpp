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
  const unsigned int degu = u.size() - 1;
  const unsigned int degv = v.size() - 1;

  // Object for product polynomial p = uv
  const unsigned int degp = degu + degv;

  Eigen::VectorXd uv(degp + 1);

  // TODO: (4-3.a) Multiply polynomials $u$ and $v$ naively.
  // START
  for (unsigned int i = 0; i <= degp; ++i) {
    const int fst = std::max(0, i - degv);
    const int lst = std::min(degu, i);
    for (int j = fst; j <= lst; ++j) {
      uv(i) += u(j) * v(i - j);
    }
  }
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
  Eigen::FFT<double> fft;
  Eigen::VectorXcd u_tmp_ = fft.fwd(u_tmp);
  Eigen::VectorXcd v_tmp_ = fft.fwd(v_tmp);
  Eigen::VectorXcd tmp = u_tmp_.cwiseProduct(v_tmp_);
  uv = fft.inv(tmp).real();
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

  // zero padding
  Eigen::VectorXd uv_tmp = uv;
  uv_tmp.conservativeResizeLike(Eigen::VectorXd::Zero(dim));
  Eigen::VectorXd u_tmp = u;
  u_tmp.conservativeResizeLike(Eigen::VectorXd::Zero(dim));

  Eigen::FFT<double> fft;
  Eigen::VectorXcd uv_tmp_ = fft.fwd(uv_tmp);
  Eigen::VectorXcd u_tmp_ = fft.fwd(u_tmp);

  // check divisibility: case(i)
  for (unsigned int i = 0; i < dim; ++i) {
    if (abs(uv_tmp_(i)) < 1e-13) {
      if (abs(u_tmp_(i)) < 1e-13) {
        // make cwiseQuotient at i-th position equal 0
        uv_tmp_(i) = 0.0;  // complex assignment (0., 0.)
        u_tmp_(i) = 1.0;   // complex assignment (1., 0.)
      } else {
        std::cerr << "uv can't be divided by u\n";
        return Eigen::VectorXd(0);
      }
    }
  }

  Eigen::VectorXcd tmp = uv_tmp_.cwiseQuotient(u_tmp_);

  v = fft.inv(tmp).real();
  // check divisibility: case(ii)
  for (unsigned int i = n; i < dim; ++i) {
    if (abs(v(i)) > 1e-13) {
      std::cerr << "uv can't be divided by u\n";
      return Eigen::VectorXd(0);
    }
  }

  // reshape v to a suitable size
  v.conservativeResizeLike(Eigen::VectorXd::Zero(n));
  // (mn-1) - (m-1) + 1
  // END

  return v;
}
/* SAM_LISTING_END_2 */
