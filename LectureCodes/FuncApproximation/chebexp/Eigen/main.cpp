# include <iostream>
# include "chebexp.hpp"

int main () {
  VectorXd y = VectorXd::LinSpaced(10,0,1);
  std::cout << "Chebychev expansion coefficients for y = [0:0.1:1]: \n" << chebexp(y) << "\n";
  return 0;
}
