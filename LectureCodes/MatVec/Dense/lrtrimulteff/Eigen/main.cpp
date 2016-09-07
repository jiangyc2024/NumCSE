#include <Eigen/Dense>
#include <figure/figure.hpp>
#include <iostream>
#include <cassert>
using namespace Eigen;
#include "lrtrimulteff.hpp"
int main(){
  int n = 3;
  MatrixXd A(n,n);
  A << 1,2,3,4,5,6,7,8,9;
  MatrixXd B(n,n);
  B << 9,8,7,6,5,4,3,2,1;
  VectorXd x(n);
  x << 4,5,6;
  VectorXd y(n), y_eff(n);
/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd AB = A*B.transpose();
y = AB.triangularView<Eigen::Upper>()*x;
/* SAM_LISTING_END_0 */
  lrtrimulteff(A, B, x, y_eff);
  std::cout << "y orginal\n" << y << std::endl;
  std::cout << "y eff\n" << y_eff << std::endl;
  return 0;
}
