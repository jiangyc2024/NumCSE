
#include "extrapolpi.hpp"
#include <math.h>

using namespace Eigen;

int main() {
  
  std::cout << std::setprecision(6);
  std::cout << "extrapolate_to_pi(5) = " << extrapolate_to_pi(5) << std::endl;
  std::cout << "extrapolate_to_pi(15) = " << extrapolate_to_pi(15) << std::endl;
  
  plotExtrapolationError();
  
  return 0;
}
