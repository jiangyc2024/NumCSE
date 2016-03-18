# include <iostream>
# include <cmath>
# include <complex>
# include <functional>
# include "./fourcoeffcomp.hpp"

int main () {
  std::complex<double> I(0,1); // imaginary unit
  std::function<std::complex<double>(double)> c = [I](double x) { 
    return std::exp(I*x*M_PI); 
  };

  std::vector<std::complex<double>> y;
  fourcoeffcomp(y, c, 4);
 
  // test the output
  for (auto v : y){
    std::cout << v << " ";
  }
  
  std::cout << "\n";
  return 0;
}
