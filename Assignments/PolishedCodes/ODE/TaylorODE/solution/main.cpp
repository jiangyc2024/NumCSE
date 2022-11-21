#include "taylorintegrator.hpp"

int main() {
  // Verify rate of convergence
  double convRate = TestCvgTaylorMethod();
  std::cout << "Estimated rate of convergence: " << convRate << std::endl;
  std::cout << "Expected rate of convergence: " << 3.0 << std::endl;
  return 0;
}
