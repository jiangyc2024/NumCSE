#ifndef NUMCSE_PCHI_HPP
#define NUMCSE_PCHI_HPP

#include <Eigen/Dense>
#include <fstream>
#include <vector>

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
/*!
 * \brief Implements a piecewise cubic Hermite interpolation.
 * Uses equidistant meshes and various methods of slope reconstruction.
 *
 */
class CubicHermiteInterpolant {
 public:
  /*!
   *! \brief Construct the slopes from the data.
   *! Use finite-differences or setting $s'(x_j) = 0$.
   *! \param[in] t Vector of nodes (assumed equidistant and sorted).
   *! \param[in] y Vector of values at nodes $t$.
   *! \param[in] s Flag to set if you want to reconstruct or set the slopes to
   *zero.
   */
  template <typename Function>
  CubicHermiteInterpolant(Function &&f, const VectorXd &t);
  virtual ~CubicHermiteInterpolant(void) = default;

  /*!
   *! \brief Evaluate the intepolant at the nodes $x$.
   *! The input is assumed sorted, unique and inside the interval
   *! \param[in] x The vector of points $x_i$ where to compute $s(x_i)$.
   *! \return Values of interpolant at $x$ (vector)
   */
  VectorXd eval(const VectorXd &x) const;

 private:
  // Provided nodes and values $(t,y)$ to compute spline,
  // Eigen vectors, $c$ contains slopes
  // All have the same size $n$
  VectorXd t_, y_, c_;
  // Difference $t(i)-t(i-1)$
  double h_;
  // Size of $t$, $y$ and $c$.
  int n_;
};
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename Function>
CubicHermiteInterpolant::CubicHermiteInterpolant(Function &&f,
                                                 const VectorXd &t)
    : t_(t), y_(t.unaryExpr(f)), c_(t.size()) {
  // Sanity check
  n_ = t_.size();
  assert(n_ == y_.size() && "t and y must have same dimension.");
  assert(n_ >= 3 && "Need at least three nodes.");
  h_ = t_(1) - t_(0);

  //// Reconstruction of the slope,
  /* SAM_LISTING_BEGIN_2 */
  // TO DO: implement reconstruction of slopes.
  // START
  // First order alternative:
  c_(0) = (-1 * y_(2) + 4 * y_(1) - 3 * y_(0)) / 2 / h_;
  for (int i = 1; i < n_ - 1; ++i) {
    c_(i) = (y_(i + 1) - y_(i - 1)) / 2 / h_;
  }

  // First order alternative:
  c_(n_ - 1) = (3 * y_(n_ - 1) - 4 * y_(n_ - 2) + 1 * y_(n_ - 3)) / 2 / h_;
  // END

  /* SAM_LISTING_END_2 */
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_3 */
VectorXd CubicHermiteInterpolant::eval(const VectorXd &x) const {
  VectorXd ret(x.size());

  // TO DO: implement evaluation function.
  // START

  int i_star = 0;
  double tmp, h_tmp, t1, t2, y1, y2, c1, c2, a1, a2, a3;
  for (int j = 0; j < x.size();) {
    // Stores the current interval index and some temporary variable
    // Find the interval for $x_j$
    // while assuming that $x$ values are sorted and within
    // $[t_{\min},t_{\max}]$
    if (t_(i_star) <= x(j) && x(j) <= t_(i_star + 1)) {
      t1 = t_(i_star);
      t2 = t_(i_star + 1);
      h_tmp = t2 - t1;  // h_tmp = h_ for equidistant grid
      y1 = y_(i_star);
      y2 = y_(i_star + 1);
      c1 = c_(i_star);
      c2 = c_(i_star + 1);
      a1 = y2 - y1;
      a2 = a1 - h_tmp * c1;
      a3 = h_tmp * c2 - a1 - a2;

      // Compute $s(x(j))$
      tmp = (x(j) - t1) / h_tmp;
      ret(j) = y1 + (a1 + (a2 + a3 * tmp) * (tmp - 1.)) * tmp;

      ++j;
    } else {
      // Otherwise, go to next interval
      ++i_star;
      // Terminate if outside any interval
      if (i_star >= t_.size() - 1) break;
      continue;
    }
  }

  // END

  return ret;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
template <typename Function>
Eigen::VectorXd fppchip(Function &&f, const Eigen::VectorXd &t,
                        const Eigen::VectorXd &x) {
  VectorXd ret{x.size()};
  // TO DO: compute ret.
  // START

  int n = x.size();
  VectorXd y = t.unaryExpr(f);
  int i_star = 0;
  double tmp, h_tmp, t1, t2, y1, y2, a1, a2, a3;
  for (int j = 0; j < x.size();) {
    // Stores the current interval index and some temporary variable
    // Find the interval for $x_j$
    // while assuming that $x$ values are sorted and within
    // $[t_{\min},t_{\max}]$
    if (t(i_star) <= x(j) && x(j) <= t(i_star + 1)) {
      t1 = t(i_star);
      t2 = t(i_star + 1);
      h_tmp = t2 - t1;  // h\_tmp = h\_ for equidistant grid
      y1 = y(i_star);
      y2 = y(i_star + 1);
      a1 = y2 - y1;
      a2 = a1;
      a3 = -a1 - a2;

      // Compute $s(x(j))$
      tmp = (x(j) - t1) / h_tmp;
      ret(j) = y1 + (a1 + (a2 + a3 * tmp) * (tmp - 1.)) * tmp;

      ++j;
    } else {
      // Otherwise, go to next interval
      ++i_star;
      // Terminate if outside any interval
      if (i_star >= t.size() - 1) break;
      continue;
    }
  }

  // END

  return ret;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
std::vector<double> fppchipConvergence(void) {
  // Interpoland
  auto f = [](double x) { return 1. / (1. + x * x); };

  double a = 5;  // Interval bounds will be (-a,a)
  int M = 1000;  // Number of  points in which to evaluate the interpoland

  // Precompute values at which evaluate f
  VectorXd x = VectorXd::LinSpaced(M, -a, a);
  VectorXd fx = x.unaryExpr(f);

  // Store error and number of nodes
  std::vector<double> N_nodes, err_zero;

  // TO DO: compute err\_zero.
  // START

  for (int i = 4; i <= 512; i = i << 1) {
    // Define subintervals and evaluate f there (find pairs (t,y))
    VectorXd t = VectorXd::LinSpaced(i, -a, a);

    // compute reconstructed slopes
    VectorXd s_zero_x = fppchip(f, t, x);

    // Compute infinity norm of error
    err_zero.push_back((s_zero_x - fx).lpNorm<Infinity>());
    N_nodes.push_back(i);

    if (i == 16) {
      plt::figure();
      plt::title("Interpolant with zero slope");
      plt::xlabel("t");
      plt::ylabel("y");
      plt::plot(x, s_zero_x, "r", {{"label", "$s_{zero}$"}});
      plt::plot(x, fx, "b--", {{"label", "$f$"}});
      plt::legend();
      plt::savefig("./cx_out/p_zero.png");
    }
  }

  // END

  return err_zero;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
std::vector<double> rspchipConververgence(void) {
  // Interpoland
  auto f = [](double x) { return 1. / (1. + x * x); };

  double a = 5;  // Interval bounds will be (-a,a)
  int M = 1000;  // Number of  points in which to evaluate the interpoland

  // Precompute values at which evaluate f
  VectorXd x = VectorXd::LinSpaced(M, -a, a);
  VectorXd fx = x.unaryExpr(f);

  // Store error and number of nodes
  std::vector<double> N_nodes, err_reconstr;

  // TO DO: compute err\_reconstr.
  // START

  for (int i = 4; i <= 512; i = i << 1) {
    // Define subintervals and evaluate f there (find pairs (t,y))
    VectorXd t = VectorXd::LinSpaced(i, -a, a);

    // Construct PCHI with zero and reconstructed slopes
    CubicHermiteInterpolant s_reconstr(f, t);

    // Evaluate interpolant
    VectorXd s_reconstr_x = s_reconstr.eval(x);

    // Compute infinity norm of error
    err_reconstr.push_back((s_reconstr_x - fx).lpNorm<Infinity>());
    N_nodes.push_back(i);

    if (i == 16) {
      plt::figure();
      plt::title("Interpolant with reconstructed slope");
      plt::xlabel("t");
      plt::ylabel("y");
      plt::plot(x, s_reconstr_x, "r", {{"label", "$s_{reconstr}$"}});
      plt::plot(x, fx, "b--", {{"label", "$f$"}});
      plt::legend();
      plt::savefig("./cx_out/p_reconstr.png");
    }
  }
  // END

  return err_reconstr;
}
/* SAM_LISTING_END_6 */

#endif  // NUMCSE_PCHI_HPP