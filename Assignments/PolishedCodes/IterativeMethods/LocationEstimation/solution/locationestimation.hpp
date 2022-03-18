/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#define _USE_MATH_DEFINES

/* SAM_LISTING_BEGIN_1 */
template <typename VECTOR, typename FUNCTION, typename JACOBIAN,
          typename RECORDER =
              std::function<void(const VECTOR &, const VECTOR &, double)>>
VECTOR gn(
    const VECTOR &init, FUNCTION &&F, JACOBIAN &&J, double rtol = 1.0E-6,
    double atol = 1.0E-8,
    RECORDER &&rec = [](const VECTOR &, const VECTOR &, double) -> void {}) {
  VECTOR x = init; // Vector for iterates $\cob{\Vx^{(k)}}$
  // Vector for Gauss-Newton correction $\cob{\Vs}$
  auto f{F(x)}; // Auxiliary vector in the image space of $\cob{F}$
  VECTOR s = J(x).householderQr().solve(f); // \Label[line]{gn:2}
  rec(x, s, f.norm());
  x = x - s; // Update of iterate
  // A posteriori termination based on absolute and relative tolerances
  while ((s.norm() > rtol * x.norm()) &&
         (s.norm() > atol)) { // \Label[line]{gn:term}
    f = F(x);
    s = J(x).householderQr().solve(f); // \Label[line]{gn:5}
    rec(x, s, f.norm());               // Possibly record iteration information
    x = x - s;                         // Update of iterate
  }
  return x;
}
/* SAM_LISTING_END_1 */

// Function for solving the location estimation problem.
// The positions of the receivers are passed in the matrix P, the arrival times
// in the vector ta. Speed of sound is set to 1
/* SAM_LISTING_BEGIN_7 */
template <typename RECORDER = std::function<
              void(const Eigen::Vector4d &, const Eigen::Vector4d &, double)>>
std::pair<Eigen::Vector3d, double> source_estimation(
    const Eigen::Matrix<double, 3, Eigen::Dynamic> &Q,
    const Eigen::VectorXd &ta,
    RECORDER &&rec = [](const Eigen::Vector4d &, const Eigen::Vector4d &,
                        double) -> void {}) {
  const int n = ta.size();
  assert(Q.cols() == n);
  // START: Student solution
  // Lambda function realizing the function $\cob{F}$ from \prbeqref{eq:F}
  auto F = [&Q, &ta](const Eigen::Vector4d &x) -> Eigen::VectorXd {
    const int n = ta.size();
    Eigen::VectorXd result(n);
    for (int j = 0; j < n; ++j) {
      result[j] = (Q.col(j) - x.head(3)).norm() - ta[j] + x[3];
    }
    return result;
  };
  // Lambda function for Jacobian $\cob{\Derv F}$ according to \prbeqref{eq:J}
  auto J = [&Q, &ta](const Eigen::Vector4d &x) -> Eigen::MatrixXd {
    const int n = ta.size();
    Eigen::MatrixXd Jac(n, 4);
    for (int j = 0; j < n; ++j) {
      const Eigen::Vector3d qp{Q.col(j) - x.head(3)};
      Jac.block(j, 0, 1, 3) = -qp.transpose() / qp.norm();
      Jac(j, 3) = 1.0;
    }
    return Jac;
  };
  // Compute initial guess for source position
  // as center of gravity of receiver locations
  Eigen::Vector3d cnt = Q.rowwise().sum() / n;
  // Initial guess t\_0 for starting time of source signal:
  Eigen::VectorXd::Index minpos;
  double min_ta = ta.minCoeff(&minpos);
  double t_0 = min_ta - (Q.col(minpos) - cnt).norm();
  
   
  // Build vector for intial guess
  Eigen::Vector4d init = (Eigen::Vector4d() << cnt, t_0).finished();
  
 
  // Gauss-Newton method with correction-based termination
  const double rtol = 1.0E-6;
  const double atol = 1.0E-8;
  Eigen::Vector4d xs = gn(init, F, J, rtol, atol, rec);
  // END: Student solution
  return {xs.head(3), xs[3]};
}
/* SAM_LISTING_END_7 */
