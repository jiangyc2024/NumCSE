#include "foursum.hpp"

#include <Eigen/Dense>
#include <iostream>


int main() {
  auto s = [](double x){ return std::sin(x)*std::cos(5*x); };
  Eigen::VectorXcd c;
  c = foursum::foursum(s, 2, 3);
  std::cout << "c:\n" << c << "\n";
  return 0;
}
