/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: October 2022
 */
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Invit {

/* SAM_LISTING_BEGIN_1 */
std::pair<double, Eigen::VectorXd> invit(Eigen::MatrixXd& A,
                                         double tol = 1.0E-6) {
  assert((A.cols() == A.rows()) && "A must be square!");
  const unsigned int n = A.cols();
  Eigen::VectorXd y_new = Eigen::VectorXd::Random(n);
  Eigen::VectorXd y_old(n);
  double l_new = y_new.transpose() * A * y_new;
  double l_old;
  const auto A_lu_dec = A.lu(); //\Label[line]{cit:1}
  y_new = y_new / y_new.norm();
  do {
    l_old = l_new;
    y_old = y_new; //\Label[line]{cit:2}
    y_new = A_lu_dec.solve(y_old); //\Label[line]{cit:3}
    y_new /= y_new.norm(); //\Label[line]{cit:4}
    l_new = y_new.transpose() * A * y_new; //\Label[line]{cit:5}
  } while (std::abs(l_new - l_old) > tol * std::abs(l_new));
  return {l_new, y_new};
}
/* SAM_LISTING_END_1 */

}  // namespace Invit

int main(int /*argc*/, char** /*argv*/) {
  Eigen::MatrixXd A = ((Eigen::Matrix3d() << 1,0.5,0.5,0.5,2,0.5,0.5,0.5,3).finished());
  auto [l,y] = Invit::invit(A);
  std::cout << "l_min = " << l << ", ev = " << y.transpose() << std::endl;
  return 0;
}
