#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include "./reshape.hpp"
int main () {
  Eigen::MatrixXd M = (10*Eigen::MatrixXd::Random(4,4)).cast<int>().cast<double>();
  reshapetest(M);
  return 0;
}
