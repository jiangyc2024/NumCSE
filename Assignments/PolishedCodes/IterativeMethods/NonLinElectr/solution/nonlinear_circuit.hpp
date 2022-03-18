#ifndef NONLINEAR_CIRCUIT_HPP
#define NONLINEAR_CIRCUIT_HPP

#include <Eigen/Dense>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//! \brief Compute the solution to the nonlinear system by using Newtons iteration
//! \param[in] alpha
//! \param[in] beta parameters
//! \param[in] Uin input voltages
//! \param[out] Uout output voltages
/* SAM_LISTING_BEGIN_0 */
void circuit(const double& alpha, const double& beta, const Eigen::VectorXd& Uin, Eigen::VectorXd& Uout) {
  constexpr double Ut = 0.5;
  constexpr double tau = 1e-6;
  const unsigned int n = Uin.size();
  
  // TODO: (8-7.b) Compute the output voltages for the given circuit
  // START
  double Uin_; // Uin of current node
  // lambda function for evaluation of F
  auto F = [alpha, beta](const Eigen::VectorXd& U, const double Uin_) {
    Eigen::VectorXd f(3);
    f << 3. * U(0) - U(1) - U(2),
         3. * U(1) - U(0) - U(2) - Uin_,
         3. * U(2) - U(0) - U(1) + alpha * (std::exp(beta * (U(2) - Uin_) / Ut) - 1);
    return f;
  };
  
  Eigen::MatrixXd J(3, 3); // the Jacobian
  J << 3, -1, -1,
      -1, 3, -1,
      -1, -1, 0; // dummy in $J(2, 2)$
  Eigen::VectorXd f(3); // the function
  
  for (unsigned int i = 0; i < n; ++i) {
    Uin_ = Uin(i);
    Eigen::VectorXd U = Eigen::VectorXd::Random(3); // random initial guess
    Eigen::VectorXd h = Eigen::VectorXd::Ones(3);
    while (h.cwiseAbs().maxCoeff() > tau * U.norm()) {
      J(2, 2) = 3. + (alpha * beta) / Ut * std::exp(beta * (U(2) - Uin_) / Ut);
      f = F(U, Uin_);
      h = J.partialPivLu().solve(f);
      U -= h;
    }
    Uout(i) = U(0);
  }
  // END
}
/* SAM_LISTING_END_0 */

//! \brief Plots the output voltage as a function of the input voltage
/* SAM_LISTING_BEGIN_1 */
void plotU() {
  constexpr unsigned int n = 20; // number of grid points
  constexpr double alpha = 8;
  constexpr double beta = 1;
  
  plt::figure();
  
  // TODO: (8-7.c) Plot Uout as a function of Uin with matplotlibcpp
  // START
  Eigen::VectorXd Uin = Eigen::VectorXd::LinSpaced(n, 0, 20);
  Eigen::VectorXd Uout(n);
  circuit(alpha, beta, Uin, Uout);
  
  plt::plot(Uin, Uout);
  plt::title("Nonlinear dependence caused by a diode");
  plt::xlabel("U_in");
  plt::ylabel("U_out");
  // END
  plt::savefig("./cx_out/Udiode.png");
}
/* SAM_LISTING_END_1 */

#endif
