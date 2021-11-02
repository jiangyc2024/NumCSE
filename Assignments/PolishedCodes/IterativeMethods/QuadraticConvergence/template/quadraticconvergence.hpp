#ifndef QUADRCONVERGENCE_HPP
#define QUADRCONVERGENCE_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

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
    
    // END
  }
  // TODO: (9-4.c) Define a member function that fetches the data contained in "info".
  // START
  
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
  
  // END
  return x;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void testSteffensen(void) {
  // TODO: (9-4.b) write a test of your implementation, that prints
  // an estimate of the zero of $f(x) = xe^x - 1$
  // START
  
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
  
  // END
}
/* SAM_LISTING_END_4 */

#endif
