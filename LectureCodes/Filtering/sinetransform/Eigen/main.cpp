#include <iostream>
#include <Eigen/Dense>
#include "sinetransform.hpp"


int main()
{
  Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(9, 0, 1);
  Eigen::VectorXd s;

  sinetransform(y, s);
  std::cout << "s:\n" << s << "\n";
  return 0;
}
