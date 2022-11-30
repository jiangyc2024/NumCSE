#include "fftreal.hpp"
#include <iostream>

int main() {
  const Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(10, 0, 1);
  Eigen::VectorXcd c;

  fftreal::fftreal(y, c);
  std::cout << "c:\n" << c << "\n";
  return 0;
}
