#include "numquad.hpp"
#include <cmath>
#include <iostream>
#include <numeric>


int main () {
  //auto f = [](double x){ return x; };
  //auto f = [](double x){ return 1/(1 + 25*x*x); };
  auto f = [](double x){ return std::sqrt(x); };
  const double a = 0;
  const double b = 1;
  //numquad(f, a, b, 20, "equidistant");
  //numquad(f, a, b, 20, "chebychev");
  std::cout << "final gauss quadrature error: 2/3 - Int_{0}^{1} sqrt(x) dx =\n"; 
  std::cout << 2./3. - numquad(f, a, b, 20, "gauss").back() << std::endl;
  std::cout << "final equidistant quadrature error: 2/3 - Int_{0}^{1} sqrt(x) dx =\n";
  std::cout << 2./3. - numquad(f, a, b, 20, "equidistant").back() << std::endl;
  std::cout << "final chebychev quadrature error: 2/3 - Int_{0}^{1} sqrt(x) dx =\n";
  std::cout << 2./3. - numquad(f, a, b, 20, "chebychev").back() << std::endl;
  return 0;
}
