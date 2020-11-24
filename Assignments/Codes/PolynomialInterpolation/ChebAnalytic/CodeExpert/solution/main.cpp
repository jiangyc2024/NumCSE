#include "chebyApprox.hpp"
#include <iomanip>
#define WIDTH 15

int main() {
  auto [err, upperbound] = compareErrorAndBound(20);
  std::cout << std::setw(6)  << std::left <<"Order"
            << std::setw(WIDTH) << std::left <<"Error"
            << std::setw(WIDTH) << std::left << "Upperbound"
            <<"\n";
  for (int i = 0; i < err.size(); i++) {
    std::cout << std::setw(6)  << std::left << i + 4
              << std::setw(WIDTH) << std::left << err[i]
              << std::setw(WIDTH) << std::left << upperbound[i]
              << "\n";
  }
  return 0;
}
