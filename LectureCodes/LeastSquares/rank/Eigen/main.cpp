#include <iostream>
#include <limits>
#include <Eigen/Dense>

using Eigen::MatrixXd;

int main() {
  MatrixXd A(3, 2);
  double eps = std::numeric_limits<double>::epsilon();
  A << 1, 1,
       sqrt(eps), 0,
       0, sqrt(eps);

  std::cout << A.fullPivLu().rank() << "\n"
            << (A.transpose() * A).fullPivLu().rank() << "\n";

  return 0;
}
