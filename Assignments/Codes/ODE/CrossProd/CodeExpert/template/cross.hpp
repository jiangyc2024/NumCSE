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
  std::vector<VectorXd> res(N + 1);
  // TO DO (13-1.e): Construct the implicit mid-point method with the class
  // implicit_RKIntegrator and execute the .solve() method. 
  // Return the vector containing all steps including initial and final value.
  // START
  // Initialize implicit RK with Butcher scheme
  
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
  std::vector<VectorXd> res(N + 1);
  // TO DO (13-1.g): Implement the linear implicit mid-point method for
  // an autonomous ODE y' = f(y), y(0) = y0. Return the vector containing 
  // all steps including initial and final value.
  // START
  // Initial step size
  
  // END
  return res;
}
/* SAM_LISTING_END_2 */




/* SAM_LISTING_BEGIN_3 */
void tab_crossprod(void) {
  // TO DO (13-1.e): solve the cross-product ODE with the implicit RK method
  // defined in solve_imp_mid. Tabulate the norms of the results at all steps.
  // START
  
  // END
  
  // TO DO (13-1.g): solve the cross-product ODE with the implicit RK method
  // defined in solve_lin_mid. Tabulate the norms of the results at all steps.
  // START
  
  // END
}
/* SAM_LISTING_END_3 */


#endif