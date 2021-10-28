#include "arccosquad.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

int main(int /*argc*/, char** /*argv*/) {

  std::cout << "\n\tStraightforward application of Gauss-Legendre quadrature"
            << std::endl;
  testConvGaussQuad();

  int n = 16;
  auto f = [] (double t) {return t*t*t;}; // I(t^3) = -5*pi/32
  std::cout << "\nf(t) = t^3: arccosWeightedQuad(f,16) = "
            << std::setprecision(16) << arccosWeightedQuad(f,n)
            << std::endl;

  auto g = [] (double t) {return 1./(t-2.);}; // I(1/(t-2)) = -1.491633131460735
  std::cout << "\ng(t) = 1/(t-2): arccosWeightedQuad(g,16) = "
            << std::setprecision(16) << arccosWeightedQuad(g,n)
            << std::endl;

  std::cout << "\n\tGauss-Legendre quadrature for transformed integral"
            << std::endl;
  testConvTrfGaussQuad();

  return 0;
}
