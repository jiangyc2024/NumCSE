# include <iostream>
# include "./recclenshaw.hpp"

int main() {
  VectorXd a(3); a << 1, 0.2, 0.5;
  double x = 0.7;

  std::cout << recclenshaw(a, x) << "\n";

  return 0;
}
