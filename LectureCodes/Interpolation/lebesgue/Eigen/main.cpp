# include <iostream>
# include "./lebesgue.hpp"

using Eigen::ArrayXd;

int main() {
  VectorXd equi = VectorXd::LinSpaced(10, -1, 1),
           cheb = ((2*ArrayXd::LinSpaced(10, 1, 10) - 1)/10*M_PI_2).cos().matrix();

  std::cout << "** Approximation of Lebesgue constant for 10 nodes in [-1,1] **\n";
  std::cout << "Equidistant: " << lebesgue(equi, 200) << "\n";
  std::cout << "Chebychev: " << lebesgue(cheb, 200) << "\n";

  return 0;
}
