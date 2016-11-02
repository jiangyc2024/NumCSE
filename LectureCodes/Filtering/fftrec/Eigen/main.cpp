# include <iostream>
# include <unsupported/Eigen/FFT>
# include "fftrec.hpp"

int main() {
  Eigen::FFT<double> fft;
  for (long n = 8; n < 1000; n *= 2) {
    VectorXcd y = VectorXcd::Random(n),c1, c2;
    c1 = fft.fwd(y);
    c2 = fftrec(y);
    std::cout << "Error at n = " << n << ": \t" 
              << (c1 - c2).norm() << "\n";
  }
  return 0;
}
