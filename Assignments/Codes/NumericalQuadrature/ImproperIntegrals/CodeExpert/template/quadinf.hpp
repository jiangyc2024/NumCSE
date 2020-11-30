#ifndef QUADINF_HPP
#define QUADINF_HPP

#include <cmath>
#include <iostream>

#include "golubwelsh.hpp"
#include "matplotlibcpp.h"

#define PI M_PI

namespace plt = matplotlibcpp;

//! @brief Compute $\int_a^b f(x) dx \approx \sum w_i f(x_i)$ (with scaling of
//! $w$ and $x$)
//! @tparam Function template type for function handle f (e.g.\ lambda function)
//! @param[in] f integrand
//! @param[in] w weights
//! @param[in] x nodes for interval $[-1,1]$
//! @param[in] a left boundary in $[a,b]$
//! @param[in] b right boundary in $[a,b]$
//! @return Approximation of integral $\int_a^b f(x) dx$

/* SAM_LISTING_BEGIN_1 */
template <class Function>
double quad(Function &&f, const Eigen::VectorXd &w, const Eigen::VectorXd &x,
            double a, double b) {
  // TO DO (8-3.e): implement transformation of quadrature rule
  // START
  return 0;
  // END
}
/* SAM_LISTING_END_1 */

//! @brief Compute $\int_{-\infty}^\infty f(x) dx$ using transformation $x =
//! \cot(t)$
//! @tparam Function template type for function handle f (e.g.\ lambda function)
//! @param[in] n number of Gauss points
//! @param[in] f integrand
//! @return Approximation of integral $\int_{-\infty}^\infty f(x) dx$

/* SAM_LISTING_BEGIN_2 */
template <class Function>
double quadinf(const int n, Function &&f) {
  // TO DO (8-3.e): define the function quadinf that integrates the input
  // lambda f over the real axis. First, trasform the integrand with a change
  // of variables and then use an n point Gauss quadrature.
  // Hint 1: lambda functions can take parameters inside the [] brackets
  // Hint 2: you may write an auxiliary function to compute the quadrature over
  //         a bounded interval.

  // START
  return 0;
  // END
}
/* SAM_LISTING_END_2 */

//! @brief perform convergence test for h(t) := exp(-(t-1)^2) and plot error

/* SAM_LISTING_BEGIN_3 */
void cvgQuadInf(void) {
  // Number of max Gauss pts.
  const int N = 100;
  plt::figure();
  // TO DO (8-3.f): plot convergence errors against number of integration nodes
  // START

  // END
  plt::savefig("./cx_out/convergence.png");
}
/* SAM_LISTING_END_3 */

#endif
