#ifndef RADIOACTIVE_HPP
#define RADIOACTIVE_HPP

#include <Eigen/Dense>
#include <vector>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * @brief Implements the Gauss-Newton method with absolute stopping criterion
 *
 * @tparam F
 * @tparam DF
 * @param[inout] x Eigen::Vector4d = (A_0, B_0, lambda_1, lambda_2), on call:
 * initial guess, on return: x s.t. f(x) = 0
 * @param f the function s.t. f(x) = 0
 * @param df the Jacobian of f
 * @param tol absolute tolerance for stopping, defaults to 1e-14
 * @return std::vector<double> the L_inf norms of the update steps
 */
/* SAM_LISTING_BEGIN_0 */
template <typename F, typename DF>
std::vector<double> GaussNewton(Eigen::Vector4d& x, F&& f, DF&& df,
                                const double tol = 1e-14) {
  std::vector<double> gn_update;

  // TODO: (8-11.d) Implement the Newton iteration for the given Least squares
  // problem. Put in gn\_update the norms of the iteration steps for further
  // investigation.
  // START
  do {
    Eigen::Vector4d s = df(x).colPivHouseholderQr().solve(f(x));
    x -= s;
    gn_update.push_back(s.lpNorm<Eigen::Infinity>());
  } while (gn_update.back() > tol);
  // END
  return gn_update;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Plots the fitted $\Phi_B$ along with the measured data. Also plots
 * the norms of the iteration updates.
 *
 * @param t Eigen::VectorXd time values
 * @param m Eigen::VectorXd corresponding measured values of $\Phi_B$
 * @param x Eigen::VectorXd s.t. F(x) = 0, i.e. the solution given by the
 * Newton iteration
 * @param gn_update std::vector<double> norms of the update steps
 */
void plot(const Eigen::ArrayXd& t, const Eigen::VectorXd& m,
          const Eigen::VectorXd& x, const std::vector<double>& gn_update) {
  plt::figure();
  // TODO: (8-11.d) Plot the fitted PhiB along with the measurements using
  // matplotlibcpp
  auto f = [&t](const Eigen::Array4d& x) -> Eigen::VectorXd {
    const double a0 = x[0], b0 = x[1], l1 = x[2], l2 = x[3];
    return (-l2 * t).exp() * b0 +
           (l1 / (l2 - l1)) * ((-l1 * t).exp() - (-l2 * t).exp()) * a0;
  };
  plt::plot(t, m, ".", {{"label", "measured PhiB"}});
  plt::plot(t, f(x), {{"label", "fitted PhiB"}});
  const Eigen::VectorXd phiA = x[0] * (-x[2] * Eigen::ArrayXd(t)).exp();
  plt::plot(t, phiA, {{"label", "fitted PhiA"}});
  plt::xlabel("t");
  plt::legend();

  plt::savefig("./cx_out/solution.png");
  plt::figure();

  // TODO: (8-11.d) Plot the iteration error using matplotlibcpp.
  Eigen::VectorXd linspace =
      Eigen::VectorXd::LinSpaced(gn_update.size(), 1, gn_update.size());
  plt::semilogy(linspace, gn_update);
  plt::grid();
  plt::xlabel("iteration");
  plt::ylabel("change");

  plt::savefig("./cx_out/convergence.png");
}

#endif
