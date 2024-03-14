/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2022
 */

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace EnMin {
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd groundState_fpit(unsigned int n, double rtol = 1.0E-8,
                                 double atol = 1.0E-10,
                                 unsigned int maxit = 20) {
  const double h = 1.0 / (1.0 + n);
  // Vector for storing the iterates
  Eigen::VectorXd u(n);
  // Initialize the matrix $\cob{\VA}$
  Eigen::MatrixXd A(n, n);
  // Abuse the vector u as temporary storage to avoid an excessive
  // number of calls to the exponential function.
  for (unsigned int j = 0; j < n; ++j) {
    u[j] = 2.0 / h * std::exp(h * j);
  }
  for (unsigned int l = 0; l < n; ++l) {
    A(l, l) = 2.0;
    for (unsigned int j = 0; j < l; ++j) {
      A(l, l) += u[l - j];
    }
    for (unsigned int j = l + 1; j < n; ++j) {
      A(l, l) += u[j - l];
    }
    for (unsigned int j = 0; j < l; ++j) {
      A(l, j) = A(j, l) = u[l - j];
    }
  }
  // \com{Precompute} LU factorization of the matrix $\cob{\VA}$, see
  // \lref{cpp:seqsolvelsesmart}
  auto A_ludec = A.lu();
  // Initial guess zero
  u.setZero();
  // Main loop for fixed-point iteration
  unsigned int it_cnt = 0;  // iteration count
  double corr_norm;         // Norm of correction
  Eigen::VectorXd b(n);     // Stores the vector $\cob{\Vb(\Vu)}$
  do {
    // Initialize vector $\cob{\Vb(\Vu)}$
    for (unsigned int j = 0; j < n; ++j) {
      b[j] = -h * std::exp(u[j]);
    }
    // Backward and forward substitution
    Eigen::VectorXd u_new = A_ludec.solve(b); // \Label[line]{enm:1}
    corr_norm = (u - u_new).norm();
    // Next iterate
    u = u_new;
  } while ((corr_norm > rtol * u.norm()) && (corr_norm > atol) &&
           (++it_cnt < maxit));
  return u;
}
/* SAM_LISTING_END_1 */

}  // namespace EnMin

int main(int /*argc*/, char** /*argv*/) {
  std::cout
      << "NumCSE code: C++ wrapper class for discrete ground state computation"
      << std::endl;

  std::cout << "Run n = 100:\n" << EnMin::groundState_fpit(100) << std::endl;

  return 0;
}
