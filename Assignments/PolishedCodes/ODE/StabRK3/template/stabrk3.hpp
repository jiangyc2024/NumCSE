#ifndef STABRK3_H_
#define STABRK3_H_

#include <Eigen/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

/* SAM_LISTING_BEGIN_0 */
Eigen::Vector2d PredPrey(Eigen::Vector2d y0, double T, unsigned int M) {
  double h = T / M;
  Eigen::Vector2d y = y0;
  // TO DO: 12-7.a
  // START

  // END
  return y;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void SimulatePredPrey() {
  // TO DO: 12-7.b
  // START
  
  // END
}
/* SAM_LISTING_END_1 */

#endif  // #ifndef STABRK3_H_
