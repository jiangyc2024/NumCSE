#ifndef CIRCLE_APPROX_HPP
#define CIRCLE_APPROX_HPP

#include <Eigen/Dense>
#include <cassert>
#include <cmath>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * @brief Computes the solution of Ax = b in least squares sense using singular
 * value decomposition
 *
 * @param A Coefficient matrix
 * @param b right hand side vector
 * @return Eigen::MatrixXd x s.t. Ax = b (in least squares sense)
 */
inline Eigen::MatrixXd lsqSVD(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& b) {
  return A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
}

/**
 * @brief Computes the solution of Ax = b in least squares sense using QR
 * decomposition (Householder)
 *
 * @param A Coefficient matrix
 * @param b right hand side vector
 * @return Eigen::MatrixXd x s.t. Ax = b (in least squares sense)
 */
inline Eigen::MatrixXd lsqHHR(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& b) {
  return A.colPivHouseholderQr().solve(b);
}

/**
 * @brief Computes the solution of Ax = b in least squares sense using the
 * normal equations and a Cholesky solver (A^T*A is s.p.d)
 *
 * @param A Coefficient matrix
 * @param b right hand side vector
 * @return Eigen::MatrixXd x s.t. Ax = b (in least squares sense)
 */
inline Eigen::MatrixXd lsqNRM(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& b) {
  return (A.transpose() * A).ldlt().solve(A.transpose() * b);
}

/**
 * @brief
 *
 * @param x
 * @param y
 * @return Eigen::Vector3d
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::Vector3d circl_alg_fit(const Eigen::VectorXd& x,
                              const Eigen::VectorXd& y) {
  Eigen::Vector3d z = Eigen::VectorXd::Zero(3);
  assert(x.size() == y.size() && "Size mismatch!");
  const unsigned int n = x.size();

  // TODO: (9-12.b) Fit the circle in the least squares sense.
  // START
  Eigen::MatrixXd A(n, 3);
  A << -2. * x, -2. * y, -Eigen::VectorXd::Ones(n);
  // use Eigen's arrays which provide componentwise operations
  Eigen::VectorXd b = -x.array() * x.array() - y.array() * y.array();

  z = lsqHHR(A, b);  // lsqSVD(A, b); lsqNRM(A, b);
  z(2) = std::sqrt(z(2) + z(0) * z(0) + z(1) * z(1));
  // END
  return z;
}
/* SAM_LISTING_END_1 */

/**
 * @brief
 *
 * @param x
 * @param y
 * @param[inout] z
 * @param[out] err
 */
/* SAM_LISTING_BEGIN_2 */
void circl_geo_fit_GN(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
                      Eigen::Vector3d& z, std::vector<double>* err = nullptr) {
  assert(x.size() == y.size() && "Size mismatch!");
  const unsigned int n = x.size();
  constexpr double tol = 1e-14;

  // TODO: (9-12.e)
  // START
  auto R = [&x, &y](const Eigen::Vector3d& z) -> Eigen::ArrayXd {
    return sqrt((x.array() - z(0)) * (x.array() - z(0)) +
                (y.array() - z(1)) * (y.array() - z(1)));
  };

  auto F = [&R](const Eigen::Vector3d& z) -> Eigen::VectorXd {
    return R(z) - z(2);
  };

  auto DF = [&x, &y, &R, &n](const Eigen::Vector3d& z) -> Eigen::MatrixX3d {
    Eigen::MatrixXd ret(n, 3);
    ret << -(x.array() - z(0)) / R, -(y.array() - z(1)) / R,
        -Eigen::VectorXd::Ones(n);
  };

  double err_norm;
  do {
    Eigen::Vector3d s =
        lsqHHR(DF(z), F(z));  // lsqSVD(DF(z), F(z)); lsqNRM(DF(z), F(z));
    z -= s;
    err_norm = s.lpNorm<Eigen::Infinity>();
    if (err) {
      err->push_back(err_norm);
    }
  } while (err_norm > tol);
  // END
}
/* SAM_LISTING_END_2 */

/**
 * @brief
 *
 * @param x
 * @param y
 * @param[inout] z
 * @param[out] err
 */
/* SAM_LISTING_BEGIN_3 */
void circl_geo_fit_N(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
                     Eigen::Vector3d& z, std::vector<double>* err = nullptr) {
  assert(x.size() == y.size() && "Size mismatch!");
  const unsigned int n = x.size();
  constexpr double tol = 1e-14;

  // TODO: (9-12.f)
  // START
  auto gradPhi = [&n](const Eigen::Vector3d& z, const Eigen::ArrayXd& xm,
                      const Eigen::ArrayXd& ym,
                      const Eigen::ArrayXd& R) -> Eigen::Vector3d {
    Eigen::Vector3d ret;
    ret << xm.matrix().dot(z(2) / R - 1), ym.matrix().dot(z(2) / R - 1),
        n * z(2) - R.sum();
    return ret;
  };

  auto HessianPhi = [&n](const Eigen::Vector3d& z, const Eigen::ArrayXd& xm,
                         const Eigen::ArrayXd& ym,
                         const Eigen::ArrayXd& R) -> Eigen::Matrix3d {
    const double sum_x = (xm / R).sum();
    const double sum_y = (ym / R).sum();
    const double mixed_sum = z(2) * (xm * ym / pow(R, 3));
    const double inv_sum = z(2) * (1 / R).sum();
    Eigen::Matrix3d ret;
    ret << n - inv_sum + z(2) * (pow(xm, 2) / pow(R, 3)).sum(), mixed_sum,
        sum_x, mixed_sum, n - inv_sum + z(2) * (pow(ym, 2) / pow(R, 3)).sum(),
        sum_y, sum_x, sum_y, n;
    return ret;
  };

  double err_norm;
  do {
    Eigen::ArrayXd xm = x.array() - z(0);
    Eigen::ArrayXd ym = y.array() - z(1);
    Eigen::ArrayXd R = sqrt(xm * xm + ym * ym);
    Eigen::Vector3d s = lsqHHR(
        HessianPhi(z, xm, ym, R),
        gradPhi(
            z, xm, ym,
            R));  // lsqSVD(HessianPhi(z, xm, ym, R), gradPhi(z, xm, ym, R));
                  // lsqNRM(HessianPhi(z, xm, ym, R), gradPhi(z, xm, ym, R));
    z -= s;
    err_norm = s.lpNorm<Eigen::Infinity>();
    if (err) {
      err->push_back(err_norm);
    }
  } while (err_norm > tol);
  // END
}
/* SAM_LISTING_END_3 */

/**
 * @brief
 *
 */
void compare_convergence(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  plt::figure();

  // TODO: (9-12.g)
  // START
  const Eigen::Vector3d z_GN = circl_alg_fit(x, y);
  const Eigen::Vector3d z_N = z_GN;

  std::vector<double> err_GN, err_N;
  circl_geo_fit_GN(x, y, z_GN, &err_GN);
  circl_geo_fit_N(x, y, z_N, &err_N);

  Eigen::VectorXd x_GN =
      Eigen::VectorXd::LinSpaced(err_GN.size(), 1, err_GN.size());
  Eigen::VectorXd x_N =
      Eigen::VectorXd::LinSpaced(err_N.size(), 1, err_n.size());

  plt::semilogy(x_GN, err_GN, {{"label", "Gauss-Newton"}});
  plt::semilogy(x_N, err_N, {{"label", "Newton"}});

  plt::legend();
  plt::grid();
  plt::xlabel("iteration");
  plt::ylabel("change");
  // END
  plt::savefig("./cx_out/convergence.eps");  //! TODO
}

/**
 * @brief
 *
 * @param x
 * @param y
 * @return Eigen::Vector3d
 */
/* SAM_LISTING_BEGIN_4 */
Eigen::Vector3d circl_svd_fit(const Eigen::VectorXd& x,
                              const Eigen::VectorXd& y) {}
/* SAM_LISTING_END_4 */

/**
 * @brief
 *
 * @param x
 * @param y
 * @param z_alg
 * @param z_geo_GN
 * @param z_geo_N
 * @param z_svd
 */
/* SAM_LISTING_BEGIN_5 */
void plot(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
          const Eigen::Vector3d& z_alg, const Eigen::Vector3d& z_geo_GN,
          const Eigen::Vector3d& z_geo_N, const Eigen::Vector3d& z_svd) {}
/* SAM_LISTING_END_5 */

#endif