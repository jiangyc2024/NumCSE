#ifndef ADAPTINTERP_HPP
#define ADAPTINTERP_HPP

#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "intpolyval.hpp"
#include "matplotlibcpp.h"


using namespace Eigen;
namespace plt = matplotlibcpp;
/*!
 * \brief adaptivepolyintp
 * \tparam Function
 * \paramin f 
 * \paramin a
 * \paramin b
 * \paramin tol
 * \paramin N
 * \param errortab, pointer to vector
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function>
VectorXd adaptivepolyintp(Function&& f, double a, double b, double tol, int N, std::vector<double> *errortab = nullptr) {
  // TO DO (7-3.a) : implement the greedy algorithm for adaptive interpolation.
  // Ignore the errortab part of this function for now.
  
  // Generate sampling points and evaluate $f$ there
  VectorXd sampling_points = VectorXd::LinSpaced(N, a, b),
                  fvals_at_sampling_points = sampling_points.unaryExpr(f);
  //START
  
  //END
  return sampling_points; // return all sampling points
}
/* SAM_LISTING_END_1 */


/* SAM_LISTING_BEGIN_2 */
void plotInterpolationError(void) {
  // Declare test functions
  auto f1 = [](double t) { return std::sin(std::exp(2*t)); };
  auto f2 = [](double t) { return std::sqrt(t)/(1 + 16*t*t); };
  // TO DO (7-3.c): generate the plots of error vs number of nodes for f1, f2
  
  //START

  //END
}
/* SAM_LISTING_END_2 */
#endif
