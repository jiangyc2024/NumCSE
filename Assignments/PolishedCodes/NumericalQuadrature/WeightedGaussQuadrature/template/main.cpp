///
/// Minimal runner for all functions to be implemented in 8-6.i,j
///
#include <cmath>

#include "quadU.hpp"

int main() {
  auto f = [](double x){ return std::pow(x, 2); }; // \in \mathcal{P}_{2}
  // quadU should give exact results (assuming exact arithmetic) for n = 3
  std::cout << "Integrating sqrt(1 - x * x) * x * x from -1 to 1: " << quadU(f, 3) << " ; should give " << M_PI / 8. << std::endl;
  
  testQuadU();
}
