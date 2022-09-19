#ifndef QUADINF_HPP
#define QUADINF_HPP

#include <cmath>

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
  double I = 0.;
  // TODO: (7-3.e) implement transformation of quadrature rule
  // START
  for (unsigned int i = 0; i < w.size(); ++i) {
    I += f((x(i) + 1) * (b - a) / 2 + a) * w(i);
  }
  I *= (b - a) / 2.;
  // END
  return I;
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
  double I = 0.;
  // TODO: (7-3.e) define the function quadinf that integrates the input
  // lambda f over the real axis. First, transform the integrand with a change
  // of variables and then use an n point Gauss quadrature.
  // Hint 1: lambda functions can take parameters inside the [] brackets
  // Hint 2: you may write an auxiliary function to compute the quadrature over
  //         a bounded interval, see function quad above.

  // START
  Eigen::VectorXd w, x;

  // Compute nodes and weights of Gauss quadrature rule
  // using Golub-Welsh algorithm
  golubwelsh(n, w, x);

  //! NOTE: no function cot available in c++, need to resort to trigonometric
  //! identities. Both lines below are valid, the first computes three
  //! trigonometric functions, but the second is prone to cancellation
  auto ftilde = [&f](double x) {
    return f(std::cos(x) / std::sin(x)) / pow(std::sin(x), 2);
  };
  /* auto ftilde = [&f] (double x) { double cot = std::tan(PI_HALF - x); return
   * f(cot) * (1. + pow(cot,2)); }; */

  // integrate ftilde from O to PI using Gauss quadrature rules
  I = quad(ftilde, w, x, 0, PI);
  // END
  return I;
}
/* SAM_LISTING_END_2 */

//! @brief perform convergence test for h(t) := exp(-(t-1)^2) and plot error

/* SAM_LISTING_BEGIN_3 */
void cvgQuadInf(void) {
  // Number of max Gauss pts.
  constexpr int N = 100;
  plt::figure();
  // TODO: (7-3.f) plot convergence errors against number of integration nodes
  // START
  // variables for the plots
  std::vector<int> num_nodes;
  std::vector<double> err;
  std::vector<double> cvg_rate;

  // Integrand and exact integral
  auto f = [](double t) { return std::exp(-std::pow((t - 1), 2)); };
  const double I_ex = std::sqrt(PI);

  for (unsigned int n = 5; n <= N; n += 5) {
    // Approximated value of integral
    const double QS = quadinf(n, f);
    // plot data
    num_nodes.push_back(n);
    err.push_back(std::abs(QS - I_ex));
    //'straight' line in plot
    cvg_rate.push_back(std::pow(0.73, n));
  }
  plt::semilogy(num_nodes, err, "+", {{"label", "Error"}});
  plt::semilogy(num_nodes, cvg_rate, "--", {{"label", "O(0.73$^n$)"}});
  plt::legend("best");
  plt::xlabel("No. of quadrature nodes");
  plt::ylabel("|Error|");

  plt::grid("True");
  // END
  plt::savefig("./cx_out/convergence.png");
}
/* SAM_LISTING_END_3 */

#endif
