/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <limits>
#include <Eigen/Dense>

using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
int main() {
  MatrixXd A(3,2);
  // Inquire about machine precision $\to$ \cref{ex:IEEEeps}
  double eps = std::numeric_limits<double>::epsilon();
  // << initialization of matrix $\to$ \cref{par:eigeninit}
  A << 1, 1, sqrt(eps), 0, 0, sqrt(eps);
  // Output rank of \Blue{$\VA^{\top}\VA$}
  std::cout << "Rank of A: " << A.fullPivLu().rank() << std::endl
            << "Rank of A^TA: "
            << (A.transpose() * A).fullPivLu().rank() << std::endl;
  return 0;
}
/* SAM_LISTING_END_0 */
