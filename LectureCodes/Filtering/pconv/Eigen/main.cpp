#include "./myconv.hpp"
#include "./pconv.hpp"
#include "./pconvfft.hpp"
#include <iostream>

int main() {
  const long n = 5;
  const VectorXcd u = VectorXcd::Random(n), x = VectorXcd::Random(n);
  std::cout << pconv(u, x).transpose() << "\n";
  std::cout << pconvfft(u, x).transpose() << "\n";
  std::cout << dconv(u, x).transpose() << "\n";
  std::cout << fastconv(u, x).transpose() << "\n";
  return 0;
}
