# include <iostream>
# include <complex>
# include "deblur.hpp"
# include "../../blur/Eigen/blur.hpp"

int main() {
  MatrixXcd X(3,3); 
  X << 1,2,3,4,5,6,7,8,9;

  std::cout << "Testing fft2 --------------------------------------\n"
            << "X:\n" << X << "\n\n"
            << "fft2(X):\n" << fft2(X) << "\n\n"
            << "ifft2(fft2(X)):\n" << ifft2(fft2(X)) << "\n\n";

  MatrixXd C = MatrixXd::Random(3,3), 
           S = MatrixXd::Random(5,5);

  std::cout << "Testing deblur with random matrices ---------------\n"
            << deblur(C, S) << "\n";

  return 0;
}
