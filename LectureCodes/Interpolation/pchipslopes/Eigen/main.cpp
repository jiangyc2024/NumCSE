# include "./pchipslopes.hpp"
# include <iostream>
# include <Eigen/Dense>

int main() {
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(10, 0, 1),
                  y = t.cwiseProduct(t),
                  c;

  pchipslopes(t, y, c);
  std::cout << c.transpose() << "\n";
  return 0;
}
