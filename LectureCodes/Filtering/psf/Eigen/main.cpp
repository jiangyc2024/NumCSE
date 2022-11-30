#include "psf.hpp"
#include <iostream>

int main() {
  Eigen::MatrixXd S;
  const long L = 2;
  psf::psf(L, S);
  std::cout << "S:\n" << S << "\n";
  return 0;
}
