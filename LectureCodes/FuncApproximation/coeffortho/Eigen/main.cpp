# include "coeffortho.hpp"
# include <Eigen/Dense>
# include <iostream>


int main() {
  Eigen::VectorXd a;
  Eigen::VectorXd b;
  Eigen::VectorXd const t = Eigen::VectorXd::LinSpaced(10, 0, 1); 
  const Eigen::Index n = 5;

  coeffortho::coeffortho(t, n, a, b);

  std::cout << a.transpose() << "\n\n" << b.transpose() << "\n";
  return 0;
}
