# include "./recclenshaw.hpp"
# include <iostream>

int main() {
  Eigen::VectorXd a(3); a << 1, 0.2, 0.5;
  const double x = 0.7;

  std::cout << recclenshaw::recclenshaw(a, x) << "\n";

  return 0;
}
