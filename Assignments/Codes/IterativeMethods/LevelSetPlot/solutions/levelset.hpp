#pragma once

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#define PI M_PI

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
template <typename Functor>
Vector2d pointLevelSet(Functor &&f, const Vector2d &d, double c,
                       const Vector2d &x0, const Vector2d &x1,
                       double rtol = 1e-10, double atol = 1e-16) {
  Vector2d intersect;
  // START
  // The function whose zero has to found 
  auto F = [&f, &d, &c](double t) -> double { return f(t * d) - c; };

  double upd = 1.;
  // initial guesses for secand method 
  double t0 = x0.norm()/d.norm();
  double t1 = x1.norm()/d.norm();

  double fo = F(t0);
  // check termination criterion 
  while (std::abs(upd) > std::max(atol, rtol * std::min(t0, t1))) {
    double fn = F(t1);
    // secant update
    upd = fn * (t1 - t0) / (fn - fo);
    // new iteration points
    t0 = t1;
    t1 = t1 - upd;
    fo = fn;
  }
  // compute intersection of the boundary of the level set with the ray
  intersect = t1 * d;
  // END
  return intersect;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename Functor>
double areaLevelSet(Functor &&f, unsigned int n, double c) {
  assert(n > 2 && "n too small");
  double area;
  
  // define initial guess
  Vector2d x0 = {1, 0};
  Vector2d x1 = {2, 0};
  // START
  // constant angle
  double sin_ct = std::sin(2 * PI / n) * 0.5;
  VectorXd side(n + 1);
  for (unsigned int j = 0; j < n; ++j) {
    Vector2d d;
    d << std::cos(2 * PI * j / n), std::sin(2 * PI * j / n);
    // compute lengths of archimedean radii
    side(j) = pointLevelSet(f, d, c, x0, x1).norm();
  }
  // periodic extension for vectorization
  side(n) = side(0);
  // multiply pairs of adjacent sides
  area = (side.tail(n).cwiseProduct(side.head(n))).sum() * sin_ct;
  // END
  return area;
}
/* SAM_LISTING_END_1 */
