# include "chebexp.hpp"
# include <iostream>

int main () {
  // Conversion of Chebychev expansion into monomail form
  Eigen::VectorXd a(16);
  for (int n=0; n < 15; ++n) {
    a = Eigen::VectorXd::Zero(16);
    a[n] = 1.0;
    std::cout << "n = " << n << " : monomial coefficients = "
	      << std::endl << chebexp::chebexpToMonom(a).transpose() << std::endl;
  }
  
  const Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(10,0,1);
  std::cout << "Chebychev expansion coefficients for y = [0:0.1:1]: \n" << chebexp::chebexp(y) << "\n";
  return 0;
}

