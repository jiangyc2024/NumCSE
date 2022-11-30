#include <iostream>

#include "LV.hpp"

int main() {
  // The test uses the input u0=2.8, v0=1.5, T=2
  std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW = PhiAndW(2.8, 1.5, 2);
  std::cout << "Test of PhiAndW():\nPhi = " << PaW.first.transpose()
            << "\nW = \n"
            << PaW.second << "\n\n";

  auto [u0, v0] = findInitCond();
  std::cout << "Test of findInitCond():\nu0 = " << u0 << "\n"
            << "v0 = " << v0 << std::endl;

  return 0;
}