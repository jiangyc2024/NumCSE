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

  // END

  /* SAM_LISTING_END_2 */
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_3 */
VectorXd CubicHermiteInterpolant::eval(const VectorXd &x) const {
  VectorXd ret(x.size());

  // TO DO: implement evaluation function.
  // START

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

  // END

  return err_reconstr;
}
/* SAM_LISTING_END_6 */

#endif  // NUMCSE_PCHI_HPP