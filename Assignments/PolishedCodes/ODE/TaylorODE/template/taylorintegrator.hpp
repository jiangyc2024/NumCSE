#ifndef TAYLORINTEGRATORHPP
#define TAYLORINTEGRATORHPP

#include <cassert>
#include <vector>

#include <iomanip>
#include <iostream>

#include <Eigen/Dense>

#include "ode45.hpp"

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Vector2d> SolvePredPreyTaylor(double alpha1, double beta1,
                                                 double alpha2, double beta2,
                                                 double T,
                                                 const Eigen::Vector2d &y0,
                                                 unsigned int M) {
  // Vector storing the states
  std::vector<Eigen::Vector2d> res;
  // TO DO: 11-8.c
  // START

  // END
  return res;
}
/* SAM_LISTING_END_1 */

void PrintErrorTable(const Eigen::ArrayXd &M, const Eigen::ArrayXd &error) {
  std::cout << std::setw(15) << "M" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;

  for (unsigned int i = 0; i < M.size(); ++i) {
    std::cout << std::setw(15) << M(i) << std::setw(15) << error(i);
    if (i > 0) {
      std::cout << std::setw(15) << std::log2(error(i - 1) / error(i));
    }
    std::cout << std::endl;
  }
}

/* SAM_LISTING_BEGIN_2 */
double TestCvgTaylorMethod() {
  double cvgrate = 0;
  // TO DO: 11-8.d
  // START

  // END
  return cvgrate;
}
/* SAM_LISTING_END_2 */

#endif
