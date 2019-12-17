#include <Eigen/Dense>

#include <iostream>
#include <iomanip>
#include <vector>
#include "polyfit.hpp"

//! \brief Solve the autonomous IVP y' = f(y), y(0) = y0 using Rosenbrock method
//! Use semi-implicit Rosenbrock method using Jacobian evaluation. Equidistant steps of size T/N.
//! \tparam Function function type for r.h.s. f
//! \tparam Jacobian function type for Jacobian df
//! \tparam StateType type of solution space y and initial data y0
//! \param[in] f r.h.s. func f
//! \param[in] df Jacobian df of f
//! \param[in] y0 initial data y(0)
//! \param[in] N number of equidistant steps
//! \param[in] T final time
//! \return vector of y_k for each step k from 0 to N
/* SAM_LISTING_BEGIN_0 */
template <class Function, class Jacobian, class StateType>
std::vector<StateType> solveRosenbrock(Function &&f, Jacobian &&df,
                                       const StateType &y0, unsigned int N, double T) {
  
  // Will contain all time steps
  std::vector<StateType> res(N+1);
  
  // TO DO: (13-2.c)
  // START
  
  res.at(0) = y0; // Push initial data
  const double h = T/N;
  const double a = 1. / (std::sqrt(2) + 2.);
  
  // Some temporary variables
  StateType k1, k2;
  Eigen::MatrixXd J, W;
  
  // Main loop: *up to N (performs N steps)*
  for(unsigned int i = 1; i <= N; ++i) {
    StateType & yprev = res.at(i-1);
    
    // Jacobian computation
    J = df(yprev);
    W = Eigen::MatrixXd::Identity(J.rows(),J.cols()) - a*h*J;
    
    // Reuse factorization for each step
    auto W_lu = W.partialPivLu();
    
    // Increments
    k1 = W_lu.solve(f(yprev));
    k2 = W_lu.solve(f(yprev+0.5*h*k1) - a*h*J*k1);
    
    // Push new step
    res.at(i) = yprev + h*k2;
  }
  
  // END
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
double cvgRosenbrock(void) {
  
  double cvgRate=0;
  // TO DO: (13-2.d)
  // START
  
  // Final time
  const double T = 10;
  // Mesh sizes h=2^{-k} for k in K.
  const Eigen::ArrayXd K = Eigen::ArrayXd::LinSpaced(7,4,10);
  
  // Initial data
  Eigen::Vector2d y0;
  y0 << 1., 1.;
  // Parameter and useful matrix for f
  const double lambda = 1;
  Eigen::Matrix2d R;
  R << 0., -1., 1., 0.;

  // Function and its Jacobian
  auto f =  [&R, &lambda] (const Eigen::Vector2d & y) { return R*y + lambda*(1. - y.squaredNorm())*y; };
  auto df = [&lambda] (const Eigen::Vector2d & y) {
    double x = 1 - y.squaredNorm();
    Eigen::Matrix2d J;
    J << lambda*x - 2*lambda*y(0)*y(0),
    -1 - 2*lambda*y(1)*y(0),
    1 - 2*lambda*y(1)*y(0),
    lambda*x - 2*lambda*y(1)*y(1);
    return J;
  };
  
  // Reference mesh size
  const int N_ref = 10*std::pow(2,12);  
  // Reference solution
  auto solref = solveRosenbrock(f, df, y0, N_ref, T);
  
  Eigen::ArrayXd Error(K.size());
  std::cout << std::setw(15) << "N" << std::setw(16) << "maxerr\n";
  // Main loop: loop over all meshes
  for(unsigned int i = 0; i < K.size(); ++i) {
    // h = 2^{-k} => N = T*h = T*2^k
    int N = T*std::pow(2,K[i]);
    // Get solution
    auto sol = solveRosenbrock(f, df, y0, N, T);
    // Compute error
    double maxerr = 0;
    for(unsigned int j = 0; j < sol.size(); ++j) {
        maxerr = std::max(maxerr, (sol.at(j) - solref.at((j*N_ref)/N)).norm());
    }
    
    Error[i] = maxerr;
    std::cout << std::setw(15) << N << std::setw(16) << maxerr << std::endl;
  }
  // Use log(N)=log(T*2^k)=log(T)+log(2)*k to get natural logarithm of N.
  cvgRate = -polyfit(std::log(2)*K,Error.log(),1)(0);
  // END
  return cvgRate;
}
/* SAM_LISTING_END_1 */
