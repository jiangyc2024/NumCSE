#include <iomanip>

#include "extrapolpi.hpp"

int main() {
  std::cout << std::setprecision(6);
  std::cout << "extrapolate_to_pi(5) = " << extrapolate_to_pi(5) << std::endl;
  std::cout << "extrapolate_to_pi(15) = " << extrapolate_to_pi(15) << std::endl;
  
  std::cout << "\nPlotting error for k = 15:" << std::endl;
  plotExtrapolationError(15);
  
  std::cout << "\nPlotting error for k = 30:" << std::endl;
  plotExtrapolationError(30);

  return 0;
}
