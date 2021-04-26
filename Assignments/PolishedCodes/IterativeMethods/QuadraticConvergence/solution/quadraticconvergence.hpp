#ifndef QUADRCONVERGENCE_HPP
#define QUADRCONVERGENCE_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/* SAM_LISTING_BEGIN_2 */
template <class T>
class Logger {
  std::vector<T> info;

public:
  // Hint: there is no need to implement a constructor
  // since info is initialized as empty vector by default.

  // TODO: (9-4.c) overload the operator() so that it adds its argument to
  // the member "info".
  void operator()(T val) {
    // START
    info.push_back(val);
    // END
  }
  // TODO: (9-4.c) Define a member function that fetches the data contained in "info".
  // START
  std::vector<T> getInfo(void) { return info; }
  // END

  void print_log() {
    for (auto v : info) {
      std::cout << v << std::endl;
    }
  }
};
/* SAM_LISTING_END_2 */

/*! @brief Steffensen's method
 *! @param[in] f Function handler
 *! @param[in] x0 Initial guess
 *! @param[out] x Final estimation returned by the Steffensen's method
 */
/* SAM_LISTING_BEGIN_0 */
template <class Function>
double steffensen(Function &&f, double x0) {
  double x_old = x0;
  double x = x0;
  // TODO: (9-4.a) implement the Steffensen's method for a function f
  // START
  double upd = 1;
  constexpr double eps = std::numeric_limits<double>::epsilon();
  // Iterate until machine precision is reached
  while (std::abs(upd) > eps * std::abs(x_old)) {
    x_old = x;              // Storing the old iterate
    const double fx = f(x); // Only 2 evaluations of $f$ at each step
    if (fx != 0) {
      upd = fx * fx / (f(x + fx) - fx);
      x -= upd; // New iterate
    } else {
      upd = 0;
    }
  }
  // END
  return x;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void testSteffensen(void) {
  // TODO: (9-4.b) write a test of your implementation, that prints
  // an estimate of the zero of $f(x) = xe^x - 1$
  // START
  const double x =
      steffensen([](double x) { return x * std::exp(x) - 1; }, 1.0);
  std::cout << "The iterative method converges to " << x << std::endl;
  // END
}
/* SAM_LISTING_END_1 */

/*! @brief Steffensen's method
 *! @param[in] f Function handler
 *! @param[in] x0 Initial guess
 *! @param[out] x Final estimation returned by the Steffensen's method
 */
/* SAM_LISTING_BEGIN_3 */
template <class Function>
double steffensen_log(Function &&f, double x0,
                      Logger<double> *logger_p = nullptr) {

  double x = x0;
  // TODO: (9-4.c) Modify the function steffensen: use the class Logger to
  // save all the iterations x of the Steffensen's method.
  // START
  bool log_enabled = false;
  if (logger_p != nullptr) {
    log_enabled = true;
    (*logger_p)(x); // alternative: logger_p->operator()(x)
  }
  double upd = 1;
  constexpr double eps = std::numeric_limits<double>::epsilon();
  // Iterate until machine precision is reached
  while (std::abs(upd) > eps * x) {
    const double fx = f(x); // Only 2 evaluations of $f$ at each step
    if (fx != 0) {
      upd = fx * fx / (f(x + fx) - fx);
      x -= upd;
      if (log_enabled) {
        //(*logger_p)(x); // alternative: logger_p->operator()(x)
        logger_p->operator()(x);
      }
    } else {
      upd = 0;
    }
  }
  // END
  return x;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
void orderSteffensen(void) {
  auto f = [](double x) { return x * std::exp(x) - 1; };
  constexpr double x_star = 0.567143290409784; // use as exact value
  // TODO: (9-4.c) tabulate values from which you can read the
  // order of Steffensen's method.
  // Hint: to approximate the convergence rate, use the formula
  // $(\log(e_i) - \log(e_{i-1}))/ (\log(e_{i-1}) - \log(e_{i-2}))$
  // START
  Logger<double> logger;

  steffensen_log(f, 1.0, &logger);
  std::vector<double> myData = logger.getInfo();
  const unsigned n = myData.size();
  
  plt::figure();
  plt::title("Quadratic convergence");
  plt::xlabel("ratios");
  plt::ylabel("iterates");
  std::vector<double> iterates(n - 2), ratios_(n - 2), two(n - 2);

  Eigen::VectorXd errs(n), log_errs(n);
  for (unsigned i = 0; i < n; ++i) {
    errs(i) = std::abs(myData[i] - x_star);
    log_errs(i) = std::log(errs(i));
  }
  Eigen::VectorXd ratios = Eigen::VectorXd::Zero(n);
  for (unsigned i = 2; i < n; ++i) {
    ratios(i) =
        (log_errs(i) - log_errs(i - 1)) / (log_errs(i - 1) - log_errs(i - 2));
    ratios_[i - 2] = ratios(i);
    iterates[i - 2] = i;
    two[i - 2] = 2;
  }
  
  plt::plot(iterates, ratios_, "+r", {{"label", "Steffensen"}});
  plt::plot(iterates, two, "-", {{"label", "rate = 2"}});
  plt::legend("best");
  plt::savefig("QuadraticConvergence.eps");
  // Print output
  std::cout.precision(10);
  std::cout << std::setw(20) << "x" << std::setw(20) << "errors"
            << std::setw(20) << "ratios" << std::endl;
  std::cout << "       -----------------------------------------------------"
            << std::endl;
  for (unsigned i = 0; i < n; ++i) {
    std::cout << std::setw(20) << myData[i] << std::setw(20) << errs(i)
              << std::setw(20) << ratios(i) << std::endl;
  }
  // END
}
/* SAM_LISTING_END_4 */

#endif
