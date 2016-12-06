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
    // Will contain all steps, reserve memory for efficiency
    std::vector<VectorXd> res;
    /* SAM_LISTING_BEGIN_2 */
#if SOLUTION
    // Initial step size
    double h = T / N;
    int d = y0.size();

    res.reserve(N+1);
    // Store initial data
    res.push_back(y0);

    // Initialize some memory to store temporary values
    VectorXd ytemp1 = y0;
    VectorXd ytemp2 = y0;
    // Pointers to swap previous value
    VectorXd * yold = &ytemp1;
    VectorXd * ynew = &ytemp2;
    MatrixXd eye = MatrixXd::Identity(3,3);

    // Loop over all fixed steps
    for(unsigned int k = 0; k < N; ++k) {
        // Compute, save and swap next step
        *ynew = *yold +
                h*(eye - h/2. * Jf(*yold)).lu().solve(f(*yold));
        res.push_back(*ynew);
        std::swap(yold, ynew);
    }
#else // TEMPLATE
    // TODO: implement linear implicit MPR
#endif // TEMPLATE
    return res;
    /* SAM_LISTING_END_2 */
}


int main(void) {

    // 1. Implicit mid-point method
    /* SAM_LISTING_BEGIN_1 */
    std::cout << "1. Implicit midpoint method"
              << std::endl << std::endl;

    unsigned int s = 1;
    MatrixXd A(s,s);
    VectorXd b(s);
    A << 1./2.;
    b << 1.;

#if SOLUTION
    // Coefficients and handle for ODE
    double c = 1.;
    Vector3d a;
    a << 1., 0., 0.;
    auto f = [&a, &c] (const Vector3d & y) -> Vector3d {
        return a.cross(y) + c*y.cross(a.cross(y));
    };

    auto Jf = [&a, &c] (const Vector3d & y) {
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

    // Initial value and final time
    Vector3d y0;
    y0 << 1., 1., 1.;
    double T = 10.;

    // Initialize implicit RK with Butcher scheme
    implicit_RKIntegrator RK(A,b);
    int N = 128;

    auto res = RK.solve(f, Jf, T, y0, N);
    for (int i=0; i<N+1; i++) {
        std::cout<< "norm(y(" << T*i/N << ")) = "
                 << res[i].norm() <<std::endl;
    }
#else // TEMPLATE
    // TODO: use implicit RK integrator to solve ode.
#endif
    /* SAM_LISTING_END_1 */

    // 2. Linear implicit mid-point method
    /* SAM_LISTING_BEGIN_3 */
    std::cout << std::endl
              << "2. Linear implicit midpoint method"
              << std::endl << std::endl;

#if SOLUTION
    res = solve_lin_mid(f, Jf, T, y0, N);
    for (int i=0; i<N+1; i++) {
        std::cout << "norm(y(" << T*i/N << ")) = "
                  << res[i].norm()
                  << std::endl;
    }
#else // TEMPLATE
    // TODO: use linear implicit RK integrator to solve ode.
#endif
    /* SAM_LISTING_END_3 */
}
