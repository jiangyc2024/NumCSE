#include "rkintegrator.hpp"
#include "polyfit.hpp"

//! \file rk3prey.hpp Solve prey/predator model with RK-SSM method

/* SAM_LISTING_BEGIN_0 */
double RK3prey() {
  double conv_rate = 0;
  // Dimension of state space
  unsigned int d = 2;
  // Exact value y(10) at final time T = 10 (approximated)
  Eigen::VectorXd yex(d);
  yex << 0.319465882659820, 9.730809352326228;
  
  // TO DO: Use an RKIntegrator to solve the predator/prey IVP up to T=10,
  // as described in (12-2.b). Tabulate the error as the number of steps
  // increases N = 2^7, 2^8, ..., 2^14, and estimate the convergence rate.
  // HINT: You may use polyfit() to calculate the convergence rate.
  // START
  
  // Implementation of butcher scheme
  unsigned int s = 3;
  Eigen::MatrixXd A(s,s);
  Eigen::VectorXd b(s);
  A << 0,      0,      0,
       1./3.,  0,      0,
       0,      2./3.,  0;
  b << 1./4.,  0,      3./4.;
  
  // Initialize RK with Butcher scheme
  RKIntegrator<Eigen::VectorXd> RK(A,b);
  
  // Coefficients and handle for prey/predator model
  double alpha1 = 3.;
  double alpha2 = 2.;
  double beta1 = 0.1;
  double beta2 = 0.1;
  auto f = [&alpha1, &alpha2, &beta1, &beta2] (const Eigen::VectorXd & y) { 
    auto temp = y;
    temp(0) *= alpha1 - beta1*y(1);
    temp(1) *= -alpha2 + beta2*y(0);
    return temp;
  };
  
  // Initial value for model
  Eigen::VectorXd y0(d);
  y0 << 100, 5;
  
  // Final time for model
  double T = 10.;
  
  // Array of number of steps (for convergence study)
  ArrayXd N(8);
  N << 128, 256, 512, 1024, 2048, 4096, 8192, 16384;
  ArrayXd Error(N.size());
  
  // Start convergence study
  std::cout  << std::setw(15) << "N" << std::setw(15) << "error" << std::endl;
  for(unsigned int i = 0; i < N.size(); ++i) {
    auto res = RK.solve(f, T, y0, N[i]);
    double err = (res.back() - yex).norm();
    std::cout << std::setw(15) << N[i] << std::setw(15) << err << std::endl;
    Error[i] = err;
  }
  
  VectorXd coeffs = polyfit(N.log(),Error.log(),1);
  conv_rate = -coeffs(0);
  // END
  return conv_rate;
}
/* SAM_LISTING_END_0 */
