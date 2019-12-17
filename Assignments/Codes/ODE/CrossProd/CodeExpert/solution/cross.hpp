#ifndef CROSS_HPP
#define CROSS_HPP

#include <vector>

#include "implicit_rkintegrator.hpp"

using namespace Eigen;

//! \tparam Function type for function implementing the rhs function.
//! Must have VectorXd operator()(VectorXd x)
//! \tparam Jacobian type for function implementing the Jacobian of f.
//! Must have MatrixXd operator()(VectorXd x)
//! \param[in] f function handle for rhs in y' = f(y), e.g. implemented
//! using lambda function.
//! \param[in] Jf function handle for Jf, e.g. implemented using lambda function
//! \param[in] T final time T
//! \param[in] y0 initial data y(0) = y0 for y' = f(y)
//! \param[in] N number of steps to perform.
/* SAM_LISTING_BEGIN_1 */
template <class Function, class Jacobian>
std::vector<VectorXd> solve_imp_mid(Function &&f,
                                    Jacobian &&Jf,
                                    double T,
                                    const VectorXd &y0,
                                    unsigned int N) {
  std::vector<VectorXd> res;
  // TO DO (13-1.e): Construct the implicit mid-point method with the class
  // implicit_RKIntegrator and execute the .solve() method. 
  // Return the vector containing all steps including initial and final value.
  // START
  // Initialize implicit RK with Butcher scheme
  unsigned int s = 1;
  MatrixXd A(s,s);
  VectorXd b(s);
  A << 1./2.;
  b << 1.;
  implicit_RKIntegrator RK(A,b);
  
  res = RK.solve(f, Jf, T, y0, N);
  // END
  return res;
}
/* SAM_LISTING_END_1 */



//! \tparam Function type for function implementing the rhs function.
//! Must have VectorXd operator()(VectorXd x)
//! \tparam Jacobian type for function implementing the Jacobian of f.
//! Must have MatrixXd operator()(VectorXd x)
//! \param[in] f function handle for rhs in y' = f(y), e.g. implemented
//! using lambda function.
//! \param[in] Jf function handle for Jf, e.g. implemented using lambda function
//! \param[in] T final time T
//! \param[in] y0 initial data y(0) = y0 for y' = f(y)
//! \param[in] N number of steps to perform.
/* SAM_LISTING_BEGIN_2 */
template <class Function, class Jacobian>
std::vector<VectorXd> solve_lin_mid(Function &&f,
                                    Jacobian &&Jf,
                                    double T,
                                    const VectorXd &y0,
                                    unsigned int N)  {
  std::vector<VectorXd> res;
  // TO DO (13-1.g): Implement the linear implicit mid-point method for
  // an autonomous ODE y' = f(y), y(0) = y0. Return the vector containing 
  // all steps including initial and final value.
  // START
  // Initial step size
  double h = T / N;
  int d = y0.size();
  // reserve memory for efficiency
  res.reserve(N + 1);
  // Store initial data
  res.push_back(y0);

  // Initialize some memory to store temporary values
  VectorXd ytemp1 = y0;
  VectorXd ytemp2 = y0;
  // Pointers to swap previous value
  VectorXd *yold = &ytemp1;
  VectorXd *ynew = &ytemp2;
  MatrixXd eye = MatrixXd::Identity(3,3);

  // Loop over all fixed steps
  for(unsigned int k = 0; k < N; ++k) {
      // Compute, save and swap next step
      *ynew = *yold +
              h*(eye - h/2. * Jf(*yold)).lu().solve(f(*yold));
      res.push_back(*ynew);
      std::swap(yold, ynew);
  }
  // END
  return res;
}
/* SAM_LISTING_END_2 */




/* SAM_LISTING_BEGIN_3 */
void tab_crossprod(void) {
  // TO DO (13-1.e): solve the cross-product ODE with the implicit RK method
  // defined in solve_imp_mid. Tabulate the norms of the results at all steps.
  // START
  double T = 10.;
  int N =128;
  // set data
  double c = 1.;
  Vector3d y0;
  y0 << 1., 1., 1.;
  Vector3d a;
  a << 1., 0., 0.;
  // define rhs
  auto f = [&a, &c] (const Vector3d &y) -> Vector3d {
      return a.cross(y) + c*y.cross(a.cross(y));
  };
  //define Jacobian of rhs
  auto Jf = [&a, &c] (const Vector3d &y) {
    Matrix3d temp;
    temp << -c*(a(1)*y(1) + a(2)*y(2)),
            c*(2*a(0)*y(1) - a(1)*y(0)) - a(2),
            a(1) + c*(2*a(0)*y(2) - a(2)*y(0)),
            a(2) - c*(a(0)*y(1) - 2*a(1)*y(0)),
            -c*(a(0)*y(0) + a(2)*y(2)),
            c*(2*a(1)*y(2) - a(2)*y(1)) - a(0),
            - a(1) - c*(a(0)*y(2) - 2*a(2)*y(0)),
            a(0) - c*(a(1)*y(2) - 2*a(2)*y(1)),
            -c*(a(0)*y(0) + a(1)*y(1));
    return temp;
  };
 
  std::vector<VectorXd> res_imp = solve_imp_mid(f, Jf, T, y0, N);
  
  std::cout << "1. Implicit midpoint method" << std::endl;
  std::cout << std::setw(10) << "t" 
            << std::setw(15) << "norm(y(t))" << std::endl;
  
  for (int i=0; i<N+1; ++i) {
    std::cout << std::setw(10) << T*i/N  
              << std::setw(15) << res_imp[i].norm() << std::endl;
  }
  // END
  
  // TO DO (13-1.g): solve the cross-product ODE with the implicit RK method
  // defined in solve_lin_mid. Tabulate the norms of the results at all steps.
  // START
  std::vector<VectorXd> res_lin = solve_lin_mid(f, Jf, T, y0, N);
  
  std::cout << "\n2. Linear implicit midpoint method" << std::endl;
  std::cout << std::setw(10) << "t" 
            << std::setw(15) << "norm(y(t))" << std::endl;
  for (int i=0; i<N+1; ++i) {
    std::cout << std::setw(10) << T*i/N  
              << std::setw(15) << res_lin[i].norm() << std::endl;
  }
  // END
}
/* SAM_LISTING_END_3 */


#endif