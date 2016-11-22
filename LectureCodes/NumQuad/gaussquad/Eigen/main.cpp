# include <iostream>
# include "gaussquad.hpp"

int main() {
  QuadRule qr;
  for (unsigned n = 2; n <= 8; n *= 2) {
    qr = gaussquad_(n);
    std::cout << "n = " << n << " -----------------------------------------------\n"
              << "x = " << qr.nodes.transpose() << "\n"
              << "w = " << qr.weights.transpose() << "\n";
  }
  return 0;
}
