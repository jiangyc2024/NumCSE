# include <iostream>
# include "freqfilter.hpp"

int main() {
  Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(10, 0, 1);
  Eigen::VectorXd low; Eigen::VectorXd high;
  int k = 2;

  freqfilter(y, k, low, high);
  std::cout << "low:\n" << low << "\n";
  std::cout << "high:\n" << high << "\n";
  return 0;
}
