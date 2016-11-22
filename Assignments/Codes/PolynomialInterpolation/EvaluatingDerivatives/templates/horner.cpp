#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include "timer.h"

# include <figure/figure.hpp>

/*
 * @brief Evaluate a polynomial and its derivative using Horner scheme
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] pair containing p(x),p'(x)
 */
/* SAM_LISTING_BEGIN_0 */
template <typename CoeffVec>
std::pair<double, double> evaldp (const CoeffVec& c, const double x) {
  std::pair<double, double> p;
  double px, dpx;
  int s = c.size();

    // TODO: evaluate a polynomial using Horner scheme

  return p;
}
/* SAM_LISTING_END_0 */

/*
 * @brief Evaluate a polynomial and its derivative using a naive implementation
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] pair containing p(x),p'(x)
 */
/* SAM_LISTING_BEGIN_1 */
template <typename CoeffVec>
std::pair<double, double> evaldp_naive(const CoeffVec& c, const double x) {
  std::pair<double, double> p;
  double px,dpx;
  int n=c.size();
    
    // TODO: evaluate a polynomial using naive implementation

  return p;
}
/* SAM_LISTING_END_1 */

int main() {
  std::vector<double> c {3, 1, 5, 7, 9};
  double x = .123;
    
/* SAM_LISTING_BEGIN_2 */
    // TODO: check implementations and compare runtimes
/* SAM_LISTING_END_2 */
}
