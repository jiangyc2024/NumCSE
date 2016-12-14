//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <vector>

#include <Eigen/Dense>

using namespace Eigen;

//! \brief Perform 2 steps of newton method applied to F and its jacobian DF
//! \tparam Func type for function F
//! \tparam Jac type for jacobian DF of F
//! \param[in] F function F, for which F(z) = 0 is needed
//! \param[in] DF Jacobian DF of the function F
//! \param[in,out] z initial guess and final approximation for F(z) = 0
template <class Func, class Jac>
void newton2steps(const Func & F, const Jac & DF, VectorXd & z) {
    // TODO: subproblem e: two newton steps
}

//! \brief Perform a single step of the MIRK scheme applied to the scalar ODE y' = f(y)
//! \tparam Func type for function f
//! \tparam Jac type for jacobian df of f
//! \param[in] f function f, as in y' = f(y)
//! \param[in] df Jacobian df of the function f
//! \param[in] y0 previous value
//! \param[in] h step-size
//! \return value y1 at next step
template <class Func, class Jac>
double MIRKstep(const Func & f, const Jac & df, double y0, double h) {
    // Coefficients of MIRK
    const double v1  = 1;
    const double v2  = 344./2025.;
    const double d21 = -164./2025.;
    const double b1  = 37./82.;
    const double b2  = 45./82.;
    
    // TODO: subproblem f: implement MIRK step
    return 0;
}

//! \brief Solve an ODE y' = f(y) using MIRK scheme on equidistant steps
//! \tparam Func type for function f
//! \tparam Jac type for jacobian df of f
//! \param[in] f function f, as in y' = f(y)
//! \param[in] df Jacobian df of the function f
//! \param[in] y0 initial value
//! \param[in] T final time
//! \param[in] N number of steps
//! \return value approximating y(T)
template <class Func, class Jac>
double MIRKsolve(const Func & f, const Jac & df, double y0, double T, unsigned int N) {
    
    // TODO: problem g: implement MIRK solver
    return 0;
}

int main(int, char**) {
    // r.h.s
    auto f = [] (double y) -> double { return 1 + y*y; };
    // Jacobian of $f$
    auto df = [] (double y) -> double { return 2*y; };
    // Initial data
    const double y0 = 0.;
    // Final time
    const double T = 1.;
    // Exact solution at t = T = 1
    const double yex = tan(1);
    
    //// PROBLEM h: TEST
    std::cout << "*** PROBLEM h:" << std::endl;
    // TODO: problem h: solve IVP y' = f(y) up to T
}
