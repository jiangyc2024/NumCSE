# include <iostream>
# include <Eigen/Dense>
# include "clenshaw.hpp"

int main() {
  Eigen::VectorXd a(4), x = Eigen::VectorXd::LinSpaced(10, 0, 1), y;
  a << 1, 0.2, 0.5, 0.6;

  y = clenshaw(a, x);
  std::cout << y.transpose() << "\n";

  return 0;
}
