#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "fftlsq.hpp"

using namespace Eigen;


int main(int argc, char** argv) {
  // testing testNormEqMatrix
  unsigned int n = 10;
  unsigned int m = 3;
  bool test = testNormEqMatrix(n, m);
  if (test) {
    std::cout << "testNormEqMatrix passed!\n\n";
  } else {
    std::cout << "testNormEqMatrix failed!\n\n";
  }

  // Test points
  constexpr unsigned int npoints = 10;
  VectorXd d(npoints);
  d << 0.987214, 1.03579, 0.997689, 0.917471, 1.00474, 0.92209, 1.03517,
      1.08863, 0.904992, 0.956089;

  // Find best polynomial  (coefficients)
  VectorXd g;
  g = find_c(d, m);
  std::cout << "testing find_c" << std::endl;
  std::cout << g << std::endl << std::endl;
  
  // fitting the kepler orbits and plotting the solution .png figure
  fitEllipse() ; 

}
