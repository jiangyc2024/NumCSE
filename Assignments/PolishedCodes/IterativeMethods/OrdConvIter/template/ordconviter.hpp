#ifndef ORDCONVITER_HPP
#define ORDCONVITER_HPP

#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <vector>

#include "matplotlibcpp.h"
#include "meshgrid.hpp"

namespace plt = matplotlibcpp;

/* SAM_LISTING_BEGIN_0 */
void kminplot() {
  constexpr double p = 1.5;
  constexpr double C = 2.;
  const double eps_max = std::pow(C, 1. / (1. - p));
  constexpr unsigned int ngp = 100;  // number of grid points

  plt::figure();

  // TODO: (8-5.c) Plot kmin using matplotlibcpp. To that end, start by creating
  // a mesh using the provided meshgrid function in meshgrid.hpp.
  // START

  // END
  plt::savefig("./cx_out/k_min.png");
}
/* SAM_LISTING_END_0 */
#endif
