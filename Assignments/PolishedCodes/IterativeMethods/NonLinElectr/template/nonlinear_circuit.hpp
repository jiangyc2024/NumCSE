#ifndef NONLINEAR_CIRCUIT_HPP
#define NONLINEAR_CIRCUIT_HPP

#include <Eigen/Dense>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//! \brief Compute the solution to the nonlinear system by using Newton's iteration
//! \param[in] alpha
//! \param[in] beta parameters
//! \param[in] Uin input voltages
//! \param[out] Uout output voltages
/* SAM_LISTING_BEGIN_1 */
void circuit(const double& alpha, const double& beta, const Eigen::VectorXd& Uin, Eigen::VectorXd& Uout) {
  constexpr double Ut = 0.5;
  constexpr double tau = 1e-6;
  const unsigned int n = Uin.size();
  
  // TODO: (9-7.b) Compute the output voltages for the given circuit
  // START
  
  // END
}
/* SAM_LISTING_END_1 */

//! \brief Plots the output voltage as a function of the input voltage
/* SAM_LISTING_BEGIN_1 */
void plotU() {
  constexpr unsigned int n = 20; // number of grid points
  constexpr double alpha = 8;
  constexpr double beta = 1;
  
  plt::figure();
  
  // TODO: (9-7.c) Plot Uout as a function of Uin with matplotlibcpp
  // START
  
  // END
  plt::savefig("./cx_out/Udiode.png");
}
/* SAM_LISTING_END */

#endif
