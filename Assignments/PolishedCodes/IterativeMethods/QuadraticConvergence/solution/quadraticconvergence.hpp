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

  double upd = 1;

  // Get machine precision
  const double eps = std::numeric_limits<double>::epsilon();
  
  // Iterate until machine precision is reached
  while (std::abs(upd) > eps * std::abs(x_old)) {
    // Storing the old iterate
    x_old = x;

    // Only 2 evaluations of $f$ at each step
    const double fx = f(x);

    if (fx != 0) {
      // Compute update
      upd = fx * fx / (f(x + fx) - fx);

      // New iterate
      x -= upd;
    } else {
      upd = 0;
    }
  }

  // END

  return x;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void testSteffensen() {
  // TODO: (8-4.b) write a test of your implementation, that prints
  // an estimate of the zero of $f(x) = xe^x - 1$
  // START
  
  auto f = [](double x) { return x * std::exp(x) - 1; };
  double x0 = 1.0;

  const double x = steffensen(f, x0);

  std::cout << "The iterative method converges to " << x << std::endl;
  
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
  
  // Add x to the log. Note: This is an idle lambda function if none is
  // provided. See the default value in the function signature.
  log(x);

  double upd = 1;
  const double eps = std::numeric_limits<double>::epsilon();

  // Iterate until machine precision is reached
  while (std::abs(upd) > eps * x) {
    // Only 2 evaluations of $f$ at each step
    const double fx = f(x);

    if (fx != 0) {
      // Compute update
      upd = fx * fx / (f(x + fx) - fx);

      // New iterate
      x -= upd;

      log(x);
    } else {
      upd = 0;
    }
  }
  
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
  
  std::vector<double> store;

  auto log = [&store](double x) -> void { store.push_back(x); };

  steffensen_log(f, 1.0, log);
  
  const unsigned n = store.size();

  Eigen::VectorXd errs(n), log_errs(n);
  for (unsigned i = 0; i < n; ++i) {
    errs(i) = std::abs(store[i] - x_star);
    log_errs(i) = std::log(errs(i));
  }
  
  Eigen::VectorXd ratios = Eigen::VectorXd::Zero(n);
  for (unsigned i = 2; i < n; ++i) {
    ratios(i) =
        (log_errs(i) - log_errs(i - 1)) / (log_errs(i - 1) - log_errs(i - 2));
  }
  
  // Print output
  std::cout.precision(10);
  std::cout << std::setw(20) << "x" << std::setw(20) << "errors"
            << std::setw(20) << "ratios" << std::endl;
  std::cout << "       -----------------------------------------------------"
            << std::endl;
  for (unsigned i = 0; i < n; ++i) {
    std::cout << std::setw(20) << store[i] << std::setw(20) << errs(i)
              << std::setw(20) << ratios(i) << std::endl;
  }
  
  // END
}
/* SAM_LISTING_END_3 */

#endif
