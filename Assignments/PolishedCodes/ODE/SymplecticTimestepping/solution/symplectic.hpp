#ifndef SYMPLECTIC_HPP
#define SYMPLECTIC_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

/**
 * \brief Implements a symplectic timestep for the harmonic oscillator.
 *
 * \param tau step size
 * \param pq_j state (in-place modification)
 */
/* SAM_LISTING_BEGIN_0 */
void sympTimestep(double tau, Eigen::Vector2d &pq_j) {
  // Coefficients of the method
  const Eigen::Vector3d a{2. / 3., -2. / 3., 1.};
  const Eigen::Vector3d b{7. / 24., 3. / 4., -1. / 24.};
  // TODO: (11-10.a) Implement a symplectic timestep on the harmonic oscillator
  // ODE.
  // START
  // one step of the method
  for (unsigned int i = 0; i < 3; ++i) {
    pq_j(0) += tau * b(i) * pq_j(1);
    pq_j(1) -= tau * a(i) * pq_j(0);
  }
  // END
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void sympTimesteppingODETest() {
  constexpr unsigned int nIter = 7;   // total number of iterations
  unsigned int m;                     // number of equidistant steps
  std::vector<double> errors(nIter);  // errors vector for all approx. sols
  Eigen::Vector2d approx_sol;
  approx_sol << 0.0, 1.0;  // initial conditions
  for (int k = 0; k < nIter; k++) {
    m = 10 * std::pow(2, k);
    // TODO: (11-10.b) Evaluate the error at the final step between the approx
    // solutions as given by the symplectic method and the exact solution
    // computed from the analytic formula.
    // START
    const double tau = 2.0 * M_PI / m;
    for (unsigned int i = 0; i < m; i++) {
      sympTimestep(tau, approx_sol);
    }
    // Computing the error in the maximum norm
    errors[k] = std::abs(std::sin(2.0 * M_PI) - approx_sol[0]) +
                std::abs(std::cos(2.0 * M_PI) - approx_sol[1]);
    // END
  }
  // Printing results
  std::cout << "Convergence of Symplectic Time Stepping Method:\n";
  std::cout << "\tM\t\terr(M)\n";
  for (unsigned int k = 0; k < nIter; k++) {
    std::cout << "\t" << 10 * std::pow(2, k) << "\t\t" << errors[k]
              << std::endl;
  }
}
/* SAM_LISTING_END_1 */

/**
 * \brief Solves the Hamiltonian ODE using symplectic timestepping
 *
 * \param p0 initial p
 * \param q0 initial q
 * \param T final time
 * \param M number of time steps
 * \return Eigen::MatrixXd with columns containing the evolving states
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd simulateHamiltonianDynamics(const Eigen::VectorXd &p0,
                                            const Eigen::VectorXd &q0, double T,
                                            unsigned int M) {
  const unsigned int n = p0.size();
  Eigen::MatrixXd PQ = Eigen::MatrixXd::Zero(2 * n, M + 1);
  // Coefficients of the method
  const Eigen::Vector3d a{2. / 3., -2. / 3., 1.};
  const Eigen::Vector3d b{7. / 24., 3. / 4., -1. / 24.};
  const double tau = T / M;
  Eigen::VectorXd pj(p0), qj(q0);
  PQ.col(0) << pj, qj;
  // TODO: (11-10.e) Solve the Hamiltonian ODE using symplectic timestepping.
  // Store the states in the columns of PQ.
  // START
  for (unsigned int j = 1; j <= M; j++) {   // integrate
    for (unsigned int k = 0; k < 3; k++) {  // one step
      // f(q) = - 4 * |q|^2 * q
      pj -= tau * b(k) * 4. * qj.squaredNorm() * qj;
      // g(p) = p
      qj += tau * a(k) * pj;
      PQ.col(j) << pj, qj;
    }
  }
  // END
  return PQ;
}
/* SAM_LISTING_END_3 */

#endif