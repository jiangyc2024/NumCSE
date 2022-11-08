#include <iostream>

#include "gaussquad.hpp"

int main() {
  QuadRule qr;
  for (unsigned int n = 2; n <= 8; n *= 2) {
    qr = gaussquad(n);
    std::cout << "n = " << n
              << " -----------------------------------------------\n"
              << "x = " << qr.nodes_.transpose() << "\n"
              << "w = " << qr.weights_.transpose() << "\n";
  }
  return 0;
}
