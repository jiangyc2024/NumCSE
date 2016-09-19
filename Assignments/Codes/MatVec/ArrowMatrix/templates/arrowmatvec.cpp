#include <iostream>
#include "timer.h"

#include <Eigen/Dense>

using namespace Eigen;

/* @brief Build an "arrow matrix"
 * Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
 * @param[in] d A n-dimensional vector
 * @param[in] a A n-dimensional vector
 * @param[in] x A n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
/* SAM_LISTING_BEGIN_1 */
void efficient_arrow_matrix_2_times_x(const VectorXd & d,
                                      const VectorXd & a,
                                      const VectorXd & x,
                                      VectorXd & y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
  int n = d.size();

}
/* SAM_LISTING_END_1 */

/* @brief Build an "arrow matrix"
 * Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
 * @param[in] d A n-dimensional vector
 * @param[in] a A n-dimensional vector
 * @param[in] x A n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
/* SAM_LISTING_BEGIN_0 */
void arrow_matrix_2_times_x(const VectorXd & d, const VectorXd & a,
                              const VectorXd & x, VectorXd & y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
  int n = d.size();

  VectorXd d_head = d.head(n-1);
  VectorXd a_head = a.head(n-1);
  MatrixXd d_diag = d_head.asDiagonal();

  MatrixXd A(n,n);

  A << d_diag,             a_head,
       a_head.transpose(), d(n-1);

  y = A*A*x;
}
/* SAM_LISTING_END_0 */


void runtime_arrow_matrix() {
  // TODO: your code here
}


int main(void) {
  VectorXd a(5);
  a << 1., 2., 3., 4., 5.;
  VectorXd d(5);
  d <<1., 3., 4., 5., 6.;
  VectorXd x(5);
  x << -5., 4., 6., -8., 5.;
  VectorXd yi, ye;

  arrow_matrix_2_times_x(a,d,x,yi);
  efficient_arrow_matrix_2_times_x(a,d,x,ye);

  double err = (yi - ye).norm();

  std::cout << "Error: " << err << std::endl;


  runtime_arrow_matrix();

  double eps = std::numeric_limits<double>::denorm_min();
  exit(err < eps);
}
