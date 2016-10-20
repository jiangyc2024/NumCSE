# include "./naivedft.hpp"
# include <iostream>
# include <unsupported/Eigen/FFT>

int main() {
  Eigen::FFT<double> fft;
  for (long n = 10; n < 1000; n *= 2) {
    VectorXcd y = VectorXcd::Random(n),
              c1, c2;
    fft.fwd(c1, y);
    naivedft(y, c2);
    std::cout << "Error at n = " << n << ": \t" 
              << (c1 - c2).norm() << "\n";
  }
  return 0;
}
