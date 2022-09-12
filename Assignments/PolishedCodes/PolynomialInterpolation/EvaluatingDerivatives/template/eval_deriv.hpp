#include <Eigen/Dense>
#include <complex>
#include <iomanip>

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

  // END

  plt::savefig("./cx_out/" + filename + ".png");
}
/* SAM_LISTING_END_6 */
