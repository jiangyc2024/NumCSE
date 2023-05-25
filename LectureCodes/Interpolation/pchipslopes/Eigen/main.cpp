# include "./pchipslopes.hpp"
# include <Eigen/Dense>
# include <iostream>


int main() {
  Eigen::VectorXd const t = Eigen::VectorXd::LinSpaced(10, 0, 1);
  Eigen::VectorXd const y = t.cwiseProduct(t);
  Eigen::VectorXd c;

  pchipslopes::pchipslopes(t, y, c);
  std::cout << c.transpose() << "\n";
  return 0;
}
