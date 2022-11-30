#include "myconv.hpp"
#include "pconv.hpp"
#include "pconvfft.hpp"
#include <Eigen/Dense>
#include <iostream>

using Eigen::Index;
using Eigen::VectorXcd;

int main() {
  const Index n = 5;
  const VectorXcd u = VectorXcd::Random(n);
  const VectorXcd x = VectorXcd::Random(n);
  std::cout << pconv(u, x).transpose() << "\n";
  std::cout << pconvfft(u, x).transpose() << "\n";
  std::cout << dconv(u, x).transpose() << "\n";
  std::cout << fastconv(u, x).transpose() << "\n";
  return 0;
}
