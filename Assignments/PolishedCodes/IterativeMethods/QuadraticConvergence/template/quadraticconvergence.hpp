#ifndef QUADRCONVERGENCE_HPP
#define QUADRCONVERGENCE_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

/*! @brief Steffensen's method
 *! @param[in] f Function handler
 *! @param[in] x0 Initial guess
 *! @return x Final estimation returned by the Steffensen's method
 */
/* SAM_LISTING_BEGIN_0 */
template <class Function>
double steffensen(Function &&f, double x0) {
  double x_old = x0;
  double x = x0;

  // TODO: (8-4.a) Implement the Steffensen's method for a function f.
  // START

  // END

  return x;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void testSteffensen() {
  // TODO: (8-4.b) write a test of your implementation, that prints
  // an estimate of the zero of $f(x) = xe^x - 1$
  // START

  // END
}
/* SAM_LISTING_END_1 */

/*! @brief Steffensen's method
 *! @param[in] f Function handler
 *! @param[in] x0 Initial guess
 *! @return x Final estimation returned by the Steffensen's method
 */
/* SAM_LISTING_BEGIN_2 */
template <class Function, typename LOGGER>
double steffensen_log(Function &&f, double x0,
    LOGGER &&log = [](double) -> void {}) {

  double x = x0;

  // TODO: (8-4.c) Modify the function steffensen: use the logger to
  // save all the iterations x of the Steffensen's method.
  // START

  // END

  return x;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void orderSteffensen() {
  auto f = [](double x) { return x * std::exp(x) - 1; };

  constexpr double x_star = 0.567143290409784; // use as exact value

  // TODO: (8-4.c) Tabulate values from which you can read the
  // order of Steffensen's method.
  // Hint: To approximate the convergence rate, use the formula
  // $(\log(e_i) - \log(e_{i-1}))/ (\log(e_{i-1}) - \log(e_{i-2}))$
  // START

  // END
}
/* SAM_LISTING_END_3 */

#endif
