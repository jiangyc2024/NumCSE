#include <Eigen/Dense>

#include "./blockOps.hpp"

int main () {
  Eigen::MatrixXd M = (10*Eigen::MatrixXd::Random(4,4)).cast<int>().cast<double>();
  blockAccess(M);
  return 0;
}
