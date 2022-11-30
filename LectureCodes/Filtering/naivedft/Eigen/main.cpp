#include "naivedft.hpp"
#include <iostream>
#include <unsupported/Eigen/FFT>

using Eigen::VectorXcd;

int main() {
  Eigen::FFT<double> fft;
  for (Eigen::Index n = 10; n < 1000; n *= 2) {
    const VectorXcd y = VectorXcd::Random(n);
    VectorXcd c1;
    VectorXcd c2;
    fft.fwd(c1, y);
    naivedft::naivedft(y, c2);
    std::cout << "Error at n = " << n << ": \t" 
              << (c1 - c2).norm() << "\n";
  }
  return 0;
}
