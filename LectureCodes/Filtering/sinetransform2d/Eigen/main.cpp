#include <iostream>
#include <Eigen/Dense>
#include "sinetransform2d.hpp"

int main()
{
  Eigen::VectorXd x1 = Eigen::VectorXd::LinSpaced(4, 1, 4);
  Eigen::VectorXd x2 = Eigen::VectorXd::LinSpaced(10, 1, 10);
  Eigen::MatrixXd X = x1*x2.transpose();
  std::cout << X << std::endl << std::endl;

  Eigen::MatrixXd Y;
  sinetransform2d(X, Y);
  std::cout << Y << std::endl << std::endl;

  return 0;
}
