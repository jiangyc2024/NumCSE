#include <Eigen/Dense>
#include <complex>

#include "matplotlibcpp.h"
#include "polyfit.hpp"
#include "polyval.hpp"
#include "timer.h"

namespace plt = matplotlibcpp;

/*!
 * @brief Evaluate a polynomial and its derivative using Horner scheme
 * @param[in] Vector c of size $n$, coefficients of the polynomial p
 * @param[in] Double x, where the polynomial has to be evaluated
 * @return Pair containing p(x), p'(x)
 */
/* SAM_LISTING_BEGIN_0 */
std::pair<double, double> evaldp(const Eigen::VectorXd &c, const double x) {
  double px, dpx;
  const int s = c.size();  // degree + 1

  // TODO (5-1.a): Implement the function evaldp().

  // START
  // Standard Horner scheme
  px = c(0);
  for (int i = 1; i < s; ++i) {
    px = x * px + c(i);
  }

  // Horner scheme for derivative
  dpx = (s - 1) * c(0);
  for (int i = 1; i < s - 1; ++i) {
    dpx = x * dpx + (s - i - 1) * c(i);
  }
  // END

  return {px, dpx};
}
/* SAM_LISTING_END_0 */

/*!
 * @brief Evaluate a polynomial and its derivative using a naive implementation
 * @param[in] Vector c of size $n$, coefficients of the polynomial p
 * @param[in] Double x, where the polynomial has to be evaluated
 * @return Pair containing p(x),p'(x)
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<double, double> evaldp_naive(const Eigen::VectorXd &c,
                                       const double x) {
  double px, dpx;
  const int n = c.size();

  // TODO (5-1.b): Implement the function evaldp_naive().

  // START
  px = c[0] * std::pow(x, n - 1);
  for (int i = 1; i < n; ++i) {
    px = px + c[i] * std::pow(x, n - i - 1);
  }

  dpx = (n - 1) * c[0] * std::pow(x, n - 2);
  for (int i = 1; i < n - 1; ++i) {
    dpx = dpx + (n - i - 1) * c[i] * std::pow(x, n - i - 2);
  }
  // END

  return {px, dpx};
}
/* SAM_LISTING_END_1 */

/*!
 * @brief Validates implementation of evaldp() and evaldp_naive() and
 * measures their runtimes.
 * @param[in] Const unsigned int d, max. exponent for polynomial degree 2^d
 * @return Boolean indicating if the test did fail or was passed
 */
/* SAM_LISTING_BEGIN_2 */
bool polyTestTime(const unsigned int d) {
  bool ret = false;  // Return value
  double x = 0.123;
  std::pair<double, double> p, p_naive;
  const int repeats = 10;

  std::cout
      << std::setw(10) << "n" << std::setw(25)
      << "Horner scheme:" << std::setw(25) << "Monomial approach:"
      << "\n"
      << " ================================================================\n";

  // TODO (5-1.d): Compare the output of evaldp() and evaldp_naive() for
  // $x = 0.123$ and tabulate their runtimes for $n = 2,4,...,2^d$.

  // START
  const double TOL = 1E-8;
  const int max_n = std::pow(2, d);
  Eigen::VectorXd c = Eigen::VectorXd::LinSpaced(max_n + 1, 1, max_n + 1);

  for (int n = 2; n <= max_n; n *= 2) {
    Timer tm_slow, tm_fast;

    for (int r = 0; r < repeats; r++) {
      tm_fast.start();
      p = evaldp(c.head(n + 1), x);
      tm_fast.stop();

      tm_slow.start();
      p_naive = evaldp_naive(c.head(n + 1), x);
      tm_slow.stop();
    }

    std::cout << std::setw(10) << n << std::setw(25) << tm_fast.mean()
              << std::setw(25) << tm_slow.mean() << "\n";

    // Check correctness
    double err_p = std::abs(p.first - p_naive.first);
    double err_dp = std::abs(p.second - p_naive.second);

    if (err_p >= TOL || err_dp >= TOL) {
      std::cout << "Discrepancy at n = " << n << std::endl;

      return false;
    }
  }

  ret = true;
  // END

  return ret;
}
/* SAM_LISTING_END_2 */

/*!
 * @brief Computes the derivative of the polynomial and evaluates it at x using
 * an algorithm in the spirit of Atiken-Neville.
 * @param[in] Vector t, time
 * @param[in] Vector y, data values
 * @param[in] Vector x, polynomial arguments
 * @return Vector containing the derivatives [p'(x1), ..., p'(xm)] of the
 * polynomial at positions x.
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd dipoleval(const Eigen::VectorXd &t, const Eigen::VectorXd &y,
                          const Eigen::VectorXd &x) {
  assert(t.size() == y.size() && "t and y must have same size!");
  Eigen::VectorXd ret;
  ret.resizeLike(x);  // Return value

  // TODO (5-1.e): Use the Aitken-Neville algorithm to evaluate at x
  // the derivative of the interpolating polynomial of the data (t,y).

  // START
  // Non-vectorized version: loop over all evaluation points
  for (int i = 0; i < x.size(); ++i) {
    Eigen::VectorXd p(y);
    Eigen::VectorXd dP = Eigen::VectorXd::Zero(y.size());

    // Aitken-Neville outer loop
    for (int im = 1; im < y.size(); ++im) {
      for (int i0 = im - 1; i0 >= 0; --i0) {
        dP(i0) = (p(i0 + 1) + (x(i) - t(i0)) * dP(i0 + 1) - p(i0) -
                  (x(i) - t(im)) * dP(i0)) /
                 (t(im) - t(i0));

        p(i0) = ((x(i) - t(i0)) * p(i0 + 1) - (x(i) - t(im)) * p(i0)) /
                (t(im) - t(i0));
      }
    }

    ret(i) = dP(0);
  }
  // END

  return ret;
}
/* SAM_LISTING_END_3 */

/*!
 * @brief Computes the derivative of the polynomial and evaluates it at x using
 * an alternative, less efficient algorithm. Used to check the correctness of
 * the implementation of dipoleval().
 * @param[in] Vector t, time
 * @param[in] Vector y, data values
 * @param[in] Vector x, polynomial arguments
 * @return Vector containing the derivatives [p'(x1), ..., p'(xm)] of the
 * polynomial at positions x.
 */
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd dipoleval_alt(const Eigen::VectorXd &t,
                              const Eigen::VectorXd &y,
                              const Eigen::VectorXd &x) {
  assert(t.size() == y.size() && "t and y must have same size!");
  Eigen::VectorXd ret;  // return value
  const int n = y.size();

  // TODO (5-1.f): Use polyfit() and polyval() to evaluate at x
  // the derivative of the interpolating polynomial of the data (t,y).

  // START
  Eigen::VectorXd P = polyfit(t, y, n - 1).head(n - 1);

  for (int i = 0; i < P.size(); ++i) {
    P(i) *= n - 1 - i;
  }

  ret = polyval(P, x);
  // END

  return ret;
}
/* SAM_LISTING_END_4 */

/*!
 * @brief Testing the implementation of dipoleval()
 * @return Boolean indicating if the test was passed or failed.
 */
/* SAM_LISTING_BEGIN_5 */
bool testDipolEval() {
  bool ret = false;
  Eigen::VectorXd dPx, dPx_alt;

  // TODO (5-1.f): Compare dipoleval() and dipoleval_alt() for the data
  // given on the problem sheet.

  // START
  double TOL = 1E-8;
  const int n = 10;  // polynomial degree
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(100, -1, 4);
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n + 1, 0, 3);
  Eigen::VectorXd y = sin(t.array()).matrix();

  dPx = dipoleval(t, y, x);
  dPx_alt = dipoleval_alt(t, y, x);

  double err = (dPx - dPx_alt).norm();

  if (err < TOL) {
    ret = true;
  }
  // END

  return ret;
}
/* SAM_LISTING_END_5 */

/*!
 * @brief Plots the polynomial interpolant.
 * @param[in] Filename, string storing the name of the plot.
 */
/* SAM_LISTING_BEGIN_6 */
void plotPolyInterpolant(const std::string &filename) {
  plt::figure();
  // TODO (5-1.g): Plot the derivative of the polynomial that
  // interpolates the data from (5-1.f).
  // Hint: matplotlibcpp (here = plt) implements Python's
  // matplotlib functions such as figure(), plot(), xlabel(), title(),...
  // You can use main.cpp of the LinearDataFit problem as a reference.

  // START
  Eigen::VectorXd dPx;
  const int n = 10;  // polynomial degree
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(100, 0, 3);
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n + 1, 0, 3);
  Eigen::VectorXd y = sin(t.array()).matrix();

  dPx = dipoleval(t, y, x);

  plt::plot(x, dPx);
  plt::title("Derivative of interpolating polynomial");
  plt::xlabel("t");
  plt::ylabel("p'(t)");
  // END

  plt::savefig("./cx_out/" + filename + ".png");
}
/* SAM_LISTING_END_6 */
