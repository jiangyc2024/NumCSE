#include "chebyApprox.hpp"
#include <iomanip>
#define WIDTH 15

int main() {
  std::cout << "Testing bestBound()\n";
  auto [sol_M, sol_rho] = bestBound();
  std::cout << "M = " << sol_M << std::endl;
  std::cout << "rho = " << sol_rho << std::endl;

  auto [err, upperbound] = compareErrorAndBound(20);
  std::cout << std::setw(6) << std::left << "Order" << std::setw(WIDTH)
            << std::left << "Error" << std::setw(WIDTH) << std::left
            << "Upperbound"
            << "\n";
  for (unsigned int i = 0; i < err.size(); i++) {
    std::cout << std::setw(6) << std::left << i + 4 << std::setw(WIDTH)
              << std::left << err[i] << std::setw(WIDTH) << std::left
              << upperbound[i] << "\n";
  }
  return 0;
}
