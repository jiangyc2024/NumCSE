/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: October 2022
 */

#include <cmath>
#include <iostream>

namespace Arsinh {
/* SAM_LISTING_BEGIN_1 */
double arsinh(double x) {
  if (x > 0) {
    return std::log(std::sqrt(1 + x * x) + x);
  } else {
    return -std::log(std::sqrt(1 + x * x) - x);
  }
}
/* SAM_LISTING_END_1 */

}  // namespace Arsinh

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "Implementation of arsinh" << std::endl;
  double x;
  x = 1.5;
  std::cout << "arsinh(x) = " << Arsinh::arsinh(x)
            << " <-> std::asinh(x) = " << std::asinh(x) << std::endl;
  x = -3.5;
  std::cout << "arsinh(x) = " << Arsinh::arsinh(x)
            << " <-> std::asinh(x) = " << std::asinh(x) << std::endl;
  return 0;
}
