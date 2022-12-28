#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "hermloceval.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/* SAM_LISTING_BEGIN_0 */
std::vector<Eigen::Vector2d> closedPolygonalInterpolant(
    std::vector<Eigen::Vector2d> &Sigma, const Eigen::VectorXd &x) {
  // Note: we assume that x is sorted
  assert(x.minCoeff() >= 0 && x.maxCoeff() <= 1 && "x out of range");

  std::vector<Eigen::Vector2d> res(x.size());
  unsigned int n = Sigma.size();
  // extend Sigma with periodic condition
  Sigma.push_back(Sigma[0]);

  // define variables delta_1, ..., delta_n and lambda_0,...,lambda_n
  std::vector<Eigen::Vector2d> d(n);
  std::vector<double> delta(n);
  std::vector<double> lambda(n + 1);

  lambda[0] = 0.0;
  for (unsigned int i = 0; i < n; ++i) {
    d[i] = Sigma[i + 1] - Sigma[i];
    delta[i] = (d[i]).norm();
    lambda[i + 1] = lambda[i] + delta[i];
  }

  // resume original Sigma value
  Sigma.pop_back();

  // TODO: (5-10.a)
  // START

  // END
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Vector2d> closedHermiteInterpolant(
    std::vector<Eigen::Vector2d> &Sigma, const Eigen::VectorXd &x) {
  // Note: we assume that x is sorted
  assert(x.minCoeff() >= 0 && x.maxCoeff() <= 1 && "x out of range");
  std::vector<Eigen::Vector2d> res(x.size());

  unsigned int n = Sigma.size();
  // extend Sigma with periodic condition
  Sigma.push_back(Sigma[0]);

  // define variables delta_1,...,delta_n and lambda_0,...lambda_n
  std::vector<Eigen::Vector2d> d(n);
  std::vector<double> delta(n);
  std::vector<double> lambda(n + 1);

  lambda[0] = 0.0;
  for (unsigned int i = 0; i < n; ++i) {
    d[i] = Sigma[i + 1] - Sigma[i];
    delta[i] = d[i].norm();
    lambda[i + 1] = lambda[i] + delta[i];
  }

  // define the slopes as in the definition of the interpolating
  // closed cubic Hermite curve
  std::vector<Eigen::Vector2d> slopes(n + 1);
  // TODO: (5-10.b)
  // START

  // END

  // resume original Sigma value
  Sigma.pop_back();

  return res;
}
/* SAM_LISTING_END_1 */

// clang-format off
/* SAM_LISTING_BEGIN_2 */
template <typename CurveFunctor>
std::pair<std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector2d>> adaptedHermiteInterpolant(
  CurveFunctor &&c, unsigned int nmin,
  const Eigen::VectorXd &x, double tol = 1.0e-3) {
  // clang-format on
  // Note: we assume that x is sorted
  assert(x.minCoeff() >= 0 && x.maxCoeff() <= 1 && "x out of range");
  // will contain the result
  std::vector<Eigen::Vector2d> eval;
  std::vector<Eigen::Vector2d> Sigma;

  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(nmin + 1, 0, 1);
  // TODO: (5-10.c)
  // START

  // END
  return std::make_pair(eval, Sigma);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void plotKite(const char *filename, double tol = 1.0e-3) {
  plt::figure();
  // TODO: (5-10.d)
  // START

  // END
  plt::savefig(filename);
}
/* SAM_LISTING_END_3 */
