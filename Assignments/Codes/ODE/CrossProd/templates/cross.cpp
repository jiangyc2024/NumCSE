#include <vector>

#include "implicit_rkintegrator.hpp"

using namespace Eigen;

//! \file cross.cpp Solution for CrossProd,
//! solving cross-product ODE with implicit RK method

//! \brief Solve the ODE with the linear implicit mid-point method.
//! Solve an autonomous ODE y' = f(y), y(0) = y0, using the linear implicit
//! mid-point method. Performs N equidistant steps upto time T with initial
//! data y0.
//! \tparam Function type for function implementing the rhs function.
//! Must have VectorXd operator()(VectorXd x)
//! \tparam Function2 type for function implementing the Jacobian of f.
//! Must have MatrixXd operator()(VectorXd x)
//! \param[in] f function handle for rhs in y' = f(y), e.g. implemented
//! using lambda funciton.
//! \param[in] Jf function handle for Jf, e.g. implemented using lambda funciton
//! \param[in] T final time T
//! \param[in] y0 initial data y(0) = y0 for y' = f(y)
//! \param[in] N number of steps to perform.
//! Step size is h = T / N. Steps are equidistant.
//! \return Vector containing all steps $y^n$ (for each n) including
//! initial and final value.
template <class Function, class Jacobian>
std::vector<VectorXd> solve_lin_mid(const Function &f,
                                    const Jacobian &Jf,
                                    double T,
                                    const VectorXd & y0,
                                    unsigned int N)  {
    /* SAM_LISTING_BEGIN_2 */
    // TODO: implement linear implicit MPR
    return 0;
    /* SAM_LISTING_END_2 */
}


int main(void) {

    // 1. Implicit mid-point method
    /* SAM_LISTING_END_1 */
    std::cout << "1. Implicit midpoint method"
              << std::endl << std::endl;

    unsigned int s = 1;
    MatrixXd A(s,s);
    VectorXd b(s);
    A << 1./2.;
    b << 1.;

    // TODO: use implicit RK integrator to solve ode.
    /* SAM_LISTING_END_1 */

    // 2. Linear implicit mid-point method
    /* SAM_LISTING_END_3 */
    std::cout << std::endl
              << "2. Linear implicit midpoint method"
              << std::endl << std::endl;

    // TODO: use linear implicit RK integrator to solve ode.
    /* SAM_LISTING_END_3 */
}
