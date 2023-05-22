# include "./lebesgue.hpp"
# include <iostream>

int main() {
  const Eigen::VectorXd equi = Eigen::VectorXd::LinSpaced(10, -1, 1);
  const Eigen::VectorXd cheb = ((2*Eigen::ArrayXd::LinSpaced(10, 1, 10) - 1)/10*M_PI_2).cos().matrix();

  std::cout << "** Approximation of Lebesgue constant for 10 nodes in [-1,1] **\n";
  std::cout << "Equidistant: " << lebesgue::lebesgue(equi, 200) << "\n";
  std::cout << "Chebychev: " << lebesgue::lebesgue(cheb, 200) << "\n";

  return 0;
}
