#include <iomanip>
#include "chebyApprox.hpp"
#define WIDTH 15

int main() {
  std::cout << "Testing bestBound()\n";
  auto m_rho = bestBound();
  double m = m_rho.first;
  double rho = m_rho.second;
  std::cout << "M = " << m << std::endl;
  std::cout << "rho = " << rho << std::endl;

  auto err_up = compareErrorAndBound(20);
  std::vector<double> err = err_up.first;
  std::vector<double> upperbound = err_up.second;
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
