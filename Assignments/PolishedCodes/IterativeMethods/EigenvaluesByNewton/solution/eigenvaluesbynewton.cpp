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

namespace EigByNewton {
/* SAM_LISTING_BEGIN_1 */
std::pair<double, Eigen::VectorXd> eignewton(const Eigen::MatrixXd& A,
                                             double rtol = 1.0E-6,
                                             double atol = 1.0E-8,
                                             unsigned int maxit = 20) {
  const unsigned int n = A.cols();
  assert(A.rows() == n);
  if ((A - A.transpose()).norm() > 1.0E-16 * A.norm()) {
    throw std::runtime_error("Matrix A must be symmetric");
  }
  // Vector storing the iterates, $\cob{\left(\Vz\right)_{1:n}}$ $\leftrightarrow$ $\cob{\Vx}$, $\cob{\left(\Vz\right)_{n+1}}$ $\leftrightarrow$ $\cob{\lambda}$
  Eigen::VectorXd z = Eigen::VectorXd::Zero(n + 1);
  // Determine initial guess from maximal modulus diagonal element
  Eigen::Index max_idx;
  z[n] = A.diagonal().maxCoeff(&max_idx);
  z[max_idx] = std::sqrt(2);
  // Variables for use within the main Newton loop
  Eigen::VectorXd Fz(n + 1);         // For storing $\cob{F(\Vz^{(k)})}$
  Eigen::VectorXd s(n + 1);          // Newton correction
  Eigen::MatrixXd DF(n + 1, n + 1);  // Jacobi matrix
  double s_norm;                     // Norm of Newton correction
  unsigned int it_cnt = 0;           // Iteration counter
  // Main Newton iteration loop
  do {
    // Assemble the Jacobi matrix $\cob{\Derv F(\Vz^{(k)})}$ from its blocks
    DF.block(0, 0, n, n) = A - z[n] * Eigen::MatrixXd::Identity(n, n);
    DF.block(0, n, n, 1) = -z.head(n);
    DF.block(n, 0, 1, n) = -z.head(n).transpose();
    DF(n, n) = 0.0;
    // Evaluate $\cob{F(\Vz^{(k)})}$
    Fz.head(n) = A * z.head(n) - z[n] * z.head(n);
    Fz[n] = 1 - 0.5 * z.head(n).squaredNorm();
    // Compute Newton correction $\cob{\Vs:=\Derv F(\Vz^{(k)})^{-1}F(\Vz^{(k)})}$ by Gaussian elimination
    s = DF.lu().solve(Fz);
    s_norm = s.norm();
    // Next iterate $\cob{\Vz^(k+1)}$
    z -= s;
  } while ((s_norm > rtol * z.norm()) && (s_norm > atol) && (++it_cnt < maxit));
  return {z[n], z.head(n)};
}
/* SAM_LISTING_END_1 */

}  // namespace EigByNewton

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "NumCSE code: C++ wrapper class for eigenvalue computation by "
               "Newton's method"
            << std::endl;
  Eigen::MatrixXd A(3, 3);
  A << 5, 2, 1, 2, 4, 1, 1, 1, 3;
  std::cout << "Matrix A = \n" << A << std::endl;
  auto [lambda, x] = EigByNewton::eignewton(A);
  std::cout << "eigenvalue = " << lambda
            << ", eigenvector = " << (std::sqrt(0.5) * x.transpose())
            << std::endl;
  return 0;
}
