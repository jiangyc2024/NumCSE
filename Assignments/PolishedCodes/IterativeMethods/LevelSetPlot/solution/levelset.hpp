#ifndef LEVELSET_HPP
#define LEVELSET_HPP

#include <Eigen/Dense>
#include <cassert>
#include <cmath>

/**
 * @brief Computes the intersection point between a ray and the boundary of a
 * level set.
 *
 * @tparam Functor A C++ functor object, e.g. a lambda function. Takes an
 * Eigen::Vector2d, returns a double.
 * @param f the function under consideration
 * @param d intersecting ray
 * @param c scalar defining L_c together with f
 * @param x0 the initial guess
 * @param rtol relative tolerance for stopping the iterative method
 * @param atol absolute tolerance for stopping the iterative method
 * @return Eigen::Vector2d the intersection point of ray and boundary of level
 * set
 */
/* SAM_LISTING_BEGIN_0 */
template <typename Functor>
Eigen::Vector2d pointLevelSet(Functor&& f, const Eigen::Vector2d& d, double c,
                              const Eigen::Vector2d& x0, double rtol = 1e-10,
                              double atol = 1e-16) {
  Eigen::Vector2d intersect = Eigen::VectorXd::Zero(2);

  // TODO: (8-13.c) Use the secant method to find the intersection.
  // START
  // the function whose zero has to be found
  auto F = [&f, &d, &c](double t) { return f(t * d) - c; };

  double upd = 1.;
  // initial guesses for the secant method
  double t0 = x0.norm() / d.norm();
  double t1 = 1.1 * t0;

  double fo = F(t0);
  // check correction based termination criterion
  while (std::abs(upd) > std::max(atol, rtol * std::min(t0, t1))) {
    const double fn = F(t1);
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

/**
 * @brief Computes the "Archimedean approximation" of the area of the level set.
 *
 * @tparam Functor A C++ functor object, e.g. a lambda function. Takes an
 * Eigen::Vector2d, returns a double.
 * @param f the function under consideration
 * @param n number of discretization points
 * @param c scalar defining L_c together with f
 * @return double the approximated area of the level set
 */
/* SAM_LISTING_BEGIN_1 */
template <typename Functor>
double areaLevelSet(Functor&& f, unsigned int n, double c) {
  assert(n > 2 && "n too small");
  double area = 0.;

  // define initial guesses
  Eigen::Vector2d x0 = {1, 0};

  // TODO: (8-13.d) Compute the "Archimedean approximation" of the area of the
  // level set.
  // START
  // constant angle
  const double sin_ct = std::sin(2. * M_PI / n) * 0.5;
  Eigen::VectorXd side(n + 1);
  for (unsigned int j = 0; j < n; ++j) {
    Eigen::Vector2d d;
    d << std::cos(2. * M_PI * j / n), std::sin(2. * M_PI * j / n);
    // compute length of archimedean radii
    side(j) = pointLevelSet(f, d, c, x0).norm();
  }
  // periodic extension for vectorization
  side(n) = side(0);
  // multiply pairs of adjacent sides
  area = (side.tail(n).cwiseProduct(side.head(n))).sum() * sin_ct;
  // END
  return area;
}
/* SAM_LISTING_END_1 */

#endif