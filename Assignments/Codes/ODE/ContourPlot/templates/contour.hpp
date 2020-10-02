#include "ode45.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <limits>
#include <vector>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
template <typename GradientFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic>
computeIsolinePoints(GradientFunctor &&gradF, Eigen::Vector2d y0, double T) {
  Matrix<double, 2, Dynamic> States;
  // To do: (0-2.b)
  // START
  
  // END
  return States;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEgg() {
  // crookedEggCurve will need to be reshaped to 2*(N+1).
  Matrix<double, 2, Dynamic> crookedEggCurve;
  // To do: (0-2.c)
  // START

  // END
  return crookedEggCurve;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename FFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic>
computeIsolinePointsDQ(FFunctor &&F, Eigen::Vector2d y0, double T) {
  // States will need to be reshaped to 2*(N+1).
  Matrix<double, 2, Dynamic> States;
  // To do: (0-2.d)
  // START

  // END
  return States;
}
/* SAM_LISTING_END_2 */
