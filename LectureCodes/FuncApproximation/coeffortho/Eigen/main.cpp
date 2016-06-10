# include <iostream>
# include <Eigen/Dense>
# include "coeffortho.hpp"

int main() {
  Eigen::VectorXd a, b, t = Eigen::VectorXd::LinSpaced(10, 0, 1); 
  const unsigned n = 5;

  coeffortho(t, n, a, b);

  std::cout << a.transpose() << "\n\n" << b.transpose() << "\n";
  return 0;
}
