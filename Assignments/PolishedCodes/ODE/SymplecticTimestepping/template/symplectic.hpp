#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
void sympTimestep(double tau, Eigen::Vector2d &pq_j) {
  // TO DO: 11-10.a
  // Coefficients of the method
  const Eigen::Vector3d a{2. / 3., -2. / 3., 1.};
  const Eigen::Vector3d b{7. / 24., 3. / 4., -1. / 24.};
  // one step of the method
  // START

  // END
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void sympTimesteppingODETest() {
  // TO DO: 11-10.b
  int nIter = 7;  // total number of iterations
  unsigned int m; // number of equidistant steps
  // Evaluating the error at the final step between the approx solutions as
  // given by the symplectic method and the exact solution computed from
  // the anlytic formula.
  double errors[nIter]; // errors vector for all approx. sols
  for (int k = 0; k < nIter; k++) {
    m = 10 * std::pow(2, k);
    // START

    // END
  }
  // Printing results
  std::cout << "Convergence of Symplectic Time Stepping Method:\n";
  std::cout << "\tM\t\terr(M)\n";
  for (int k = 0; k < nIter; k++) {
    std::cout << "\t" << 10 * std::pow(2, k) << "\t\t" << errors[k]
              << std::endl;
  }
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd simulateHamiltonianDynamics(const Eigen::VectorXd &p0,
                                            const Eigen::VectorXd &q0, double T,
                                            unsigned int M) {
  int n = p0.size();
  MatrixXd PQ(2 * n, M + 1);
  // TO DO: 11-10.e
  // Coefficients of the method
  Eigen::VectorXd a(3);
  a << 2. / 3., -2. / 3., 1.;
  Eigen::VectorXd b(3);
  b << 7. / 24., 3. / 4., -1. / 24.;

  double tau = T / M;
  Eigen::VectorXd pj(p0), qj(q0);
  PQ.col(0) << pj, qj;
  // START
  
  // END
  return PQ;
}
/* SAM_LISTING_END_3 */
