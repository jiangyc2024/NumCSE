#include "fftrec.hpp"
#include <iostream>
#include <unsupported/Eigen/FFT>

using Eigen::VectorXcd;

int main() {
  Eigen::FFT<double> fft;
  for (Eigen::Index n = 8; n < 1000; n *= 2) {
    const VectorXcd y = VectorXcd::Random(n);
    VectorXcd c1;
    VectorXcd c2;
    c1 = fft.fwd(y);
    c2 = fftrec::fftrec(y);
    std::cout << "Error at n = " << n << ": \t" 
              << (c1 - c2).norm() << "\n";
  }
  return 0;
}
