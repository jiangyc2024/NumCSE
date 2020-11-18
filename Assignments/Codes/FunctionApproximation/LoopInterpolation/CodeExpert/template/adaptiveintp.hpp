#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "hermloceval.hpp"
#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;

/* SAM_LISTING_BEGIN_0 */
std::vector<Vector2d> closedPolygonalInterpolant(std::vector<Vector2d> Sigma,
                                                 const VectorXd &x) {
  // Note: we assume that x is sorted
  assert(x.minCoeff() >= 0 && x.maxCoeff() <= 1 && "x out of range");

  std::vector<Vector2d> res(x.size());
  unsigned int n = Sigma.size();
  // extend Sigma with periodic condition
  Sigma.push_back(Sigma[0]);

  // define variables delta_1, ..., delta_n and lambda_0,...,lambda_n
  std::vector<Vector2d> d(n);
  std::vector<double> delta(n);
  std::vector<double> lambda(n + 1);

  lambda[0] = 0.0;
  for (unsigned int i = 0; i < n; ++i) {
    d[i] = Sigma[i + 1] - Sigma[i];
    delta[i] = (d[i]).norm();
    lambda[i + 1] = lambda[i] + delta[i];
  }
  // TO DO (6-10.a)
  // START
  unsigned int j = 0;
  // loop over sides of polygon
  for (unsigned int i = 1; i <= n; ++i) {
    // loop over evaluation points x
    while (j < x.size() && x(j) * lambda[n] <= lambda[i]) {
      // linear interpolant formula
      double ratio = (x(j) * lambda[n] - lambda[i - 1]) / delta[i - 1];
      Vector2d coord = (1 - ratio) * Sigma[i - 1] + ratio * Sigma[i];

      res[j] = coord;
      ++j;
    }
    if (j == x.size()) break;
  }
  // END
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::vector<Vector2d> closedHermiteInterpolant(std::vector<Vector2d> Sigma,
                                               const VectorXd &x) {
  // Note: we assume that x is sorted
  assert(x.minCoeff() >= 0 && x.maxCoeff() <= 1 && "x out of range");
  std::vector<Vector2d> res(x.size());

  unsigned int n = Sigma.size();
  // extend Sigma with periodic condition
  Sigma.push_back(Sigma[0]);

  // define variables delta_1,...,delta_n and lambda_0,...lambda_n
  std::vector<Vector2d> d(n);
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
  std::vector<Vector2d> slopes(n + 1);
  // TO DO (6-10.b)
  // START
  Vector2d s =
      (delta[n - 1] / delta[0] * d[0] + delta[0] / delta[n - 1] * d[n - 1]) /
      (delta[n - 1] + delta[0]);
  slopes[0] = s / s.norm();
  for (unsigned int i = 1; i < n; ++i) {
    s = (delta[i] / delta[i - 1] * d[i - 1] + +delta[i - 1] / delta[i] * d[i]) /
        (delta[i] + delta[i - 1]);
    slopes[i] = s / s.norm();
  }
  slopes[n] = slopes[0];

  VectorXd x_scl = x * lambda[n];
  unsigned int j = 0;
  // loop over sides of polygon
  for (unsigned int i = 1; i <= n; ++i) {
    // loop over evaluation points x
    while (j < x.size() && x_scl(j) <= lambda[i]) {
      // Local cubic hermite interpolation
      double coordx =
          hermloceval(x_scl(j), lambda[i - 1], lambda[i], Sigma[i - 1](0),
                      Sigma[i](0), slopes[i - 1](0), slopes[i](0));
      double coordy =
          hermloceval(x_scl(j), lambda[i - 1], lambda[i], Sigma[i - 1](1),
                      Sigma[i](1), slopes[i - 1](1), slopes[i](1));
      // store the next point
      res[j] = {coordx, coordy};
      ++j;
    }
    if (j == x.size()) break;
  }
  // END
  return res;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename CurveFunctor>
std::pair<std::vector<Vector2d>, std::vector<Vector2d>> adaptedHermiteInterpolant(CurveFunctor &&c, unsigned int nmin, const VectorXd &x, double tol = 1.0e-3) {
  // Note: we assume that x is sorted
  assert(x.minCoeff() >= 0 && x.maxCoeff() <= 1 && "x out of range");
  // will contain the result
  std::vector<Vector2d> eval;
  std::vector<Vector2d> Sigma;

  VectorXd t = VectorXd::LinSpaced(nmin + 1, 0, 1);
  // TO DO (6-10.c)
  // START
  // vector of deviations
  VectorXd dev;
  do {
    // set current number of points
    unsigned int n = t.size() - 1;
    dev.resize(n);
    Sigma.resize(n + 1);
    // evaluate Sigma for arc-length points t_i
    for (unsigned int i = 0; i <= n; ++i) {
      Sigma[i] = c(t(i));
    }
    // define variables delta_1,...,delta_n and lambda_0,...lambda_n
    VectorXd lambda(n + 1);

    VectorXd midpt(n);
    lambda[0] = 0.0;
    for (unsigned int i = 0; i < n; ++i) {
      lambda(i + 1) = lambda(i) + (Sigma[i + 1] - Sigma[i]).norm();
      midpt(i) = (lambda(i + 1) + lambda(i)) * 0.5;
    }
    // compute linear and cubic interpolant
    std::vector<Vector2d> v =
        closedPolygonalInterpolant(Sigma, midpt / lambda[n]);
    std::vector<Vector2d> w =
        closedHermiteInterpolant(Sigma, midpt / lambda[n]);
    // evaluate deviations
    for (unsigned int i = 0; i < n; ++i) {
      dev(i) = (v[i] - w[i]).norm();
    }

    double alpha = dev.mean();
    // append new arc-length values
    std::vector<double> t_temp;
    int num = 0;
    for (unsigned int j = 0; j < n; ++j) {
      if (dev(j) > 0.9 * alpha) {
        t_temp.push_back((t(j) + t(j + 1)) * 0.5);
        num++;
      }
    }
    t.conservativeResize(n + 1 + t_temp.size());

    t.tail(t_temp.size()) = VectorXd::Map(t_temp.data(), t_temp.size());
    // sort vector of arc-length values
    std::sort(t.data(), t.data() + t.size());
  } while (dev.maxCoeff() > tol);

  eval = closedHermiteInterpolant(Sigma, x);
  // END
  return std::make_pair(eval, Sigma);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void plotKite(const char *filename, double tol = 1.0e-3) {
  // TODO (6-10.d)
  // START
  // preparing call to adaptedHermiteInterpolant
  auto c = [](double t) -> Vector2d {
    Vector2d ct;
    ct << std::cos(2 * M_PI * t) + 2. / 3. * std::cos(4 * M_PI * t),
        3. / 2. * std::sin(2 * M_PI * t);
    return ct;
  };
  unsigned int nmin = 6;
  unsigned int N = 1e5;
  // initializing x, dsicarding the 1 at the end
  VectorXd x = VectorXd::LinSpaced(N + 1, 0, 1).head(N);

  std::pair<std::vector<Vector2d>, std::vector<Vector2d>> vec =
      adaptedHermiteInterpolant(c, nmin, x, tol);

  std::vector<Vector2d> eval = vec.first;
  std::vector<Vector2d> Sigma = vec.second;
  unsigned int M = Sigma.size();
  VectorXd x_hermite(N);
  VectorXd y_hermite(N);
  VectorXd x_poly(M);
  VectorXd y_poly(M);
  // assembling x and y of Hermite Interpolation
  for (unsigned int i = 0; i < N; i++) {
    x_hermite(i) = eval.at(i)(0);
    y_hermite(i) = eval.at(i)(1);
  }
  // assembling x and y of Polygonal Interpolation
  // We do not need to call ClosedPolygonalInterpolate
  // as matplotlibcpp will already connect the points
  // with straight lines for us.
  for (unsigned int i = 0; i < M; i++) {
    x_poly(i) = Sigma.at(i)(0);
    y_poly(i) = Sigma.at(i)(1);
  }
  // making the plot
  plt::figure();
  plt::title("Adapted Interpolation with tol=" + std::to_string(tol));
  plt::xlabel("x");
  plt::ylabel("y");
  plt::plot(x_poly, y_poly, "r", {{"label", "Polygonal Interpolant"}});
  plt::plot(x_hermite, y_hermite, {{"label", "Hermite Interpolation"}});
  plt::legend();
  plt::savefig(filename);
  // END
  return;
}
/* SAM_LISTING_END_3 */
