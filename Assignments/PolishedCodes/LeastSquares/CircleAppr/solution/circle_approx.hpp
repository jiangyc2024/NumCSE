#ifndef CIRCLE_APPROX_HPP
#define CIRCLE_APPROX_HPP

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

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
  return A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
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
 * @brief Computes the algebraic fit of a circle.
 *
 * @param x vector of x_i - data points
 * @param y vector of y_i - data points
 * @return Eigen::Vector3d (m_1, m_2, r) - circle description
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::Vector3d circl_alg_fit(const Eigen::VectorXd& x,
                              const Eigen::VectorXd& y) {
  Eigen::Vector3d z = Eigen::VectorXd::Zero(3);
  assert(x.size() == y.size() && "Size mismatch!");
  const unsigned int n = x.size();

  // TODO: (8-12.b) Fit the circle in the least squares sense.
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
 * @brief Computes the geometric fit to a circle using the Gauss-Newton method.
 *
 * @param x vector of x_i - data points
 * @param y vector of y_i - data points
 * @param[inout] z (m_1, m_2, r) – circle description; on entry: the initial
 * guess, on exit: the output
 * @param[out] err (optional) the error in every iteration
 */
/* SAM_LISTING_BEGIN_2 */
void circl_geo_fit_GN(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
                      Eigen::Vector3d& z, std::vector<double>* err = nullptr) {
  assert(x.size() == y.size() && "Size mismatch!");
  const unsigned int n = x.size();
  constexpr double tol = 1e-14;

  // TODO: (8-12.e) Use the Gauss-Newton method to compute the geometric fit.
  // START
  // the radius
  auto R = [&x, &y](const Eigen::Vector3d& z) -> Eigen::ArrayXd {
    return sqrt((x.array() - z(0)) * (x.array() - z(0)) +
                (y.array() - z(1)) * (y.array() - z(1)));
  };

  auto F = [&R](const Eigen::Vector3d& z) -> Eigen::VectorXd {
    return R(z) - z(2);
  };

  auto DF = [&x, &y, &R, &n](const Eigen::Vector3d& z) -> Eigen::MatrixX3d {
    Eigen::MatrixXd ret(n, 3);
    ret << -(x.array() - z(0)) / R(z), -(y.array() - z(1)) / R(z),
        -Eigen::VectorXd::Ones(n);
    return ret;
  };

  // Gauss-Newton method with correction-based stopping criterion
  // we keep track of the error-norms in err
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
 * @brief Computes the geometric fit to a circle using the Newton method.
 *
 * @param x vector of x_i - data points
 * @param y vector of y_i - data points
 * @param[inout] z (m_1, m_2, r) – circle description; on entry: the initial
 * guess, on exit: the output
 * @param[out] err (optional) the error in every iteration
 */
/* SAM_LISTING_BEGIN_3 */
void circl_geo_fit_N(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
                     Eigen::Vector3d& z, std::vector<double>* err = nullptr) {
  assert(x.size() == y.size() && "Size mismatch!");
  const unsigned int n = x.size();
  constexpr double tol = 1e-14;

  // TODO: (8-12.f) Use the Newton method to compute the geometric fit.
  // START
  auto gradPhi = [&n](const Eigen::Vector3d& z, const Eigen::ArrayXd& xm,
                      const Eigen::ArrayXd& ym,
                      const Eigen::ArrayXd& R) -> Eigen::Vector3d {
    Eigen::Vector3d ret;
    ret << xm.matrix().dot((z(2) / R - 1).matrix()),
        ym.matrix().dot((z(2) / R - 1).matrix()), n * z(2) - R.sum();
    return ret;
  };

  auto HessianPhi = [&n](const Eigen::Vector3d& z, const Eigen::ArrayXd& xm,
                         const Eigen::ArrayXd& ym,
                         const Eigen::ArrayXd& R) -> Eigen::Matrix3d {
    const double sum_x = (xm / R).sum();
    const double sum_y = (ym / R).sum();
    const double mixed_sum = z(2) * (xm * ym / pow(R, 3)).sum();
    const double inv_sum = z(2) * (1 / R).sum();
    Eigen::Matrix3d ret;
    ret << n - inv_sum + z(2) * (pow(xm, 2) / pow(R, 3)).sum(), mixed_sum,
        sum_x, mixed_sum, n - inv_sum + z(2) * (pow(ym, 2) / pow(R, 3)).sum(),
        sum_y, sum_x, sum_y, n;
    return ret;
  };

  // Newton method with correction-based stopping criterion
  // we keep track of the error-norms in err
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
 * @brief Plots the errors of the geometric fitting methods using Gauss-Newton
 * and Newton for comparison.
 *
 * @param x vector of x_i - data points
 * @param y vector of y_i - data points
 */
/* SAM_LISTING_BEGIN_4 */
void compare_convergence(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  plt::figure();
  // TODO: (8-12.g) Plot both errors using matplotlibcpp. Don't forget to use a
  // sensible axes scaling.
  // START

  // use as initial guess the algebraic fit
  Eigen::Vector3d z_GN = circl_alg_fit(x, y);
  Eigen::Vector3d z_N = z_GN;

  // keep track of the errors
  std::vector<double> err_GN, err_N;
  circl_geo_fit_GN(x, y, z_GN, &err_GN);
  circl_geo_fit_N(x, y, z_N, &err_N);

  Eigen::VectorXd x_GN =
      Eigen::VectorXd::LinSpaced(err_GN.size(), 1, err_GN.size());
  Eigen::VectorXd x_N =
      Eigen::VectorXd::LinSpaced(err_N.size(), 1, err_N.size());

  // plot with y-axis log-scaled
  plt::semilogy(x_GN, err_GN, {{"label", "Gauss-Newton"}});
  plt::semilogy(x_N, err_N, {{"label", "Newton"}});

  plt::legend();
  plt::grid();
  plt::xlabel("iteration");
  plt::ylabel("change");
  // END
  plt::savefig("./cx_out/convergence.png");
}
/* SAM_LISTING_END_4 */

/**
 * @brief Fits a circle using the constrained method and singular value
 * decomposition.
 *
 * @param x vector of x_i - data points
 * @param y vector of y_i - data points
 * @return Eigen::Vector3d (m_1, m_2, r) - circle description
 */
/* SAM_LISTING_BEGIN_5 */
Eigen::Vector3d circl_svd_fit(const Eigen::VectorXd& x,
                              const Eigen::VectorXd& y) {
  const unsigned int N = x.size();
  assert(N == y.size() && "Size mismatch!");
  Eigen::Vector3d ret = Eigen::VectorXd::Zero(3);

  // TODO: (8-12.i) Fit the circle using the constrained method and SVD.
  // START
  const Eigen::ArrayXd x_ = x.array();
  const Eigen::ArrayXd y_ = y.array();
  Eigen::MatrixXd A(x.size(), 4);
  A << (x_.square() + y_.square()).matrix(), x, y,
      Eigen::VectorXd::Ones(x.size());

  // we need to tell Eigen to compute V
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinV);
  Eigen::MatrixXd V = svd.matrixV();

  // formulas from (8-12.h)
  ret(0) = -V(1, 3) / (2. * V(0, 3));
  ret(1) = -V(2, 3) / (2. * V(0, 3));
  ret(2) = std::sqrt((V(1, 3) * V(1, 3) + V(2, 3) * V(2, 3)) /
                         (4. * V(0, 3) * V(0, 3)) -
                     V(3, 3) / V(0, 3));
  // END
  return ret;
}
/* SAM_LISTING_END_5 */

/**
 * @brief Plots the data points and fitted circles from the different methods.
 * Also, output the circle center points and radii.
 *
 * @param x vector of x_i - data points
 * @param y vector of y_i - data points
 * @param z_alg (m_1, m_2, r) - circle description as given by algebraic fit
 * @param z_geo_GN (m_1, m_2, r) - circle description as given by geometric fit
 * using Gauss-Newton
 * @param z_geo_N (m_1, m_2, r) - circle description as given by geometric fit
 * using Newton
 * @param z_svd (m_1, m_2, r) - circle description as given by constrained fit
 * using SVD
 */
void plot(const Eigen::VectorXd& x, const Eigen::VectorXd& y,
          const Eigen::Vector3d& z_alg, const Eigen::Vector3d& z_geo_GN,
          const Eigen::Vector3d& z_geo_N, const Eigen::Vector3d& z_svd) {
  plt::figure();
  constexpr unsigned int N = 100;  // number of sample points
  // TODO: (8-12.j) Plot the data points in x, y as well as the fitted circles
  // using matplotlibcpp and output the circle centers and radii.
  plt::plot(x, y, "ok", {{"label", "data points"}});

  // sample with 100 angles
  const Eigen::ArrayXd angle =
      Eigen::VectorXd::LinSpaced(N, 0, 2. * M_PI).array();
  // transformer for x-coordinate; r*cos(angle) + m_1
  auto x_coords = [&angle](const Eigen::Vector3d& z) {
    return z(2) * angle.cos().matrix() + z(0) * Eigen::VectorXd::Ones(N);
  };
  // transformer for y-coordinate; r*sin(angle) + m_2
  auto y_coords = [&angle](const Eigen::Vector3d& z) {
    return z(2) * angle.sin().matrix() + z(1) * Eigen::VectorXd::Ones(N);
  };

  // table header
  std::cout << "(m_x," << std::setw(10) << "m_y)" << std::setw(10) << "r"
            << std::endl;

  // transform, plot and output
  Eigen::VectorXd x_circle = x_coords(z_alg);
  Eigen::VectorXd y_circle = y_coords(z_alg);
  plt::plot(x_circle, y_circle, "b", {{"label", "algebraic fit"}});
  std::cout << std::setprecision(10) << "(" << z_alg(0) << ", " << std::setw(10)
            << z_alg(1) << ") " << std::setw(10) << z_alg(2) << std::endl;

  x_circle = x_coords(z_geo_GN);
  y_circle = y_coords(z_geo_GN);
  plt::plot(x_circle, y_circle, "r--",
            {{"label", "geometric fit – Gauss-Newton"}});
  std::cout << std::setprecision(10) << "(" << z_geo_GN(0) << ", "
            << std::setw(10) << z_geo_GN(1) << ") " << std::setw(10)
            << z_geo_GN(2) << std::endl;

  x_circle = x_coords(z_geo_N);
  y_circle = y_coords(z_geo_N);
  plt::plot(x_circle, y_circle, "m:", {{"label", "geometric fit – Newton"}});
  std::cout << std::setprecision(10) << "(" << z_geo_N(0) << ", "
            << std::setw(10) << z_geo_N(1) << ") " << std::setw(10)
            << z_geo_N(2) << std::endl;

  x_circle = x_coords(z_svd);
  y_circle = y_coords(z_svd);
  plt::plot(x_circle, y_circle, "c-.", {{"label", "constrained fit – SVD"}});
  std::cout << std::setprecision(10) << "(" << z_svd(0) << ", " << std::setw(10)
            << z_svd(1) << ") " << std::setw(10) << z_svd(2) << std::endl;

  plt::legend();
  plt::xlabel("x");
  plt::ylabel("y");

  plt::savefig("./cx_out/comparison.png");
}

#endif
