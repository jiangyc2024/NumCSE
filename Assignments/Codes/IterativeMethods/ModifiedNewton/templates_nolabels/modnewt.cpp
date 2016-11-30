//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <utility>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include <Eigen/Dense>

using namespace Eigen;

double norm(double x) { return std::abs(x); }
double norm(const VectorXd & x) { return x.norm(); }

/*!
 *! \brief Performs many steps of an iteration and terminate when convergence reached
 *! or maximum number of iterations has been reached.
 *! \tparam StepFunction type for the step function handle
 *! \tparam Vector argument type passed to the iteration function
 *! \tparam ErrorFunction type for the error function computing error of the method
 *! \param[in] step Function implementing the step x_{k+1} = step(x_{k}), signature Vector(const Vector&)
 *! \param[in,out] x initial data (as input) and final iteration (as output)
 *! \param[in] errf function implementing the norm of the error (errf(x)) for termination condition
 *! \param[in] eps tolerance to break iterations when res(x) < eps
 *! \param[in] max_itr maximal number of iterations
 */
template <class StepFunction, class Vector, class ErrorFunction>
bool general_nonlinear_solver(const StepFunction& step, 
                              Vector & x, 
                              const ErrorFunction& errf, 
                              double eps = 1e-8, int max_itr = 100) {
    // Temporary where to store new step
    Vector x_new = x;
    double r = 1;
    
    for(int itr = 0; itr < max_itr; ++itr) {
        
        // Compute error (or residual)
        r = errf(x);
        
        std::cout << "[Step " << itr << "] Error: " << r << std::endl;
        
        // Advance to next step, $x_{new}$ becomes $x_{k+1}$
        x_new = step(x);
        
        // Termination conditions
        // If tol reached, the we have convergence
        if (r < eps * norm(x)) {
            std::cout << "[CONVERGED] in " << itr << " it. due to err. err = " 
                      << r << " < " << eps << "." << std::endl;
            return true;
        }
        x = x_new;
    }
    
    // If max it reached
    std::cout << "[NOT CONVERGED] due to MAX it. = " << max_itr 
              << " reached, err = " << r << "." << std::endl;
    return false;
}

/*!
 *! \brief Implements a single step of the modified newton
 *! \tparam Scalar type of argument to function f: such as double, etc...
 *! \tparam Function type for the function f, likely a lambda function
 *! \tparam Jacobian type for the jacobian df, likely a lambda function
 *! \param[in] x previous value to use in step, also initial guess if needed
 *! \param[in] f function handle for f(x) = 0
 *! \param[in] df function handle for jacobian df of f
 *! \return next step x_{k+1} of modified Newton
 */
template <typename Scalar, class Function, class Jacobian>
Scalar mod_newt_step_scalar(const Scalar& x, 
                            const Function& f, 
                            const Jacobian& df) {
    Scalar y = x + f(x) / df(x);
    return y - f(y) / df(x);
}

/*!
 *! \brief Implements a single step of the modified newton
 *! \tparam Vector type of argument to function f: such as double or vector etc...
 *! \tparam Function type for the function f, likely a lambda function
 *! \tparam Jacobian type for the jacobian df, likely a lambda function
 *! \param[in] x previous value to use in step, also initial guess if needed
 *! \param[in] f function handle for f(x) = 0
 *! \param[in] df function handle for jacobian df of f
 *! \return x_next next step x_{k+1} of modified Newton
 */
template <typename Vector, class Function, class Jacobian>
Vector mod_newt_step_system(const Vector & x, 
                            const Function& f, const Jacobian& df) {
    auto lu = df(x).lu();
    // reusing LU decomposition
    Vector y = x + lu.solve(f(x));
    return y - lu.solve(f(y));
}

/**
 *! \brief Solve a scalar non-linear eq. with the modified Newton.
 * Test the empirical order of convergence of the method.
 */
void mod_newt_ord() {
    // Setting up values, functions and jacobian
    const double a = 0.123;
    auto f = [&a] (double x) { return atan(x) - a; }; // function, must see a
    auto df = [] (double x) { return 1. / (x*x+1.); }; // its derivative
    
    const double x_ex = tan(a); // exact solution
    
    // Store solution and error at each time step
    std::vector<double> sol;
    std::vector<double> err;
    
    // Compute error and push back to err, used in 
    // general\_nonlinear\_solver as breaking condition errf(x) < eps
    auto errf = [&err, x_ex] (double & x) { 
        double e =  std::abs(x - x_ex); 
        err.push_back(e); 
        return e;
    };
    
    // Perform convergence study with Modified newton for scalar
    std::cout << std::endl 
              << "*** Modified Newton method (scalar) ***"
              << std::endl << std::endl;
    std::cout << "Exact: " << x_ex << std::endl;
    
    // Initial guess and compute initial error
    double x_scalar = 5;
    sol.push_back(x_scalar);
    errf(x_scalar);
    
    // Lambda performing the next step, used to define a proper 
    // function handle to be passed to general\_nonlinear\_solver
    auto newt_scalar_step = [&sol, &f, &df] (double x) -> double { 
        double x_new = mod_newt_step_scalar(x, f, df);
        sol.push_back(x_new);
        return x_new;
    };
    
    // Actually perform the solution
    general_nonlinear_solver(newt_scalar_step, x_scalar, errf);
    
    // Print solution (final)
    std::cout << std::endl 
              << "x^*_scalar = " << x_scalar 
              << std::endl;

    // Print table of solutions, errors and EOOC
    auto space = std::setw(15);
    std::cout << space << "sol." 
              << space << "err." 
              << space << "order" 
              << std::endl;
    for(unsigned i = 0; i < sol.size(); ++ i) {
        std::cout << space << sol.at(i) 
                  << space << err.at(i);
        if(i >= 3) {
            std::cout << space << (log(err.at(i)) - log(err.at(i-1))) / (log(err.at(i-1)) - log(err.at(i-2)));
        }
        std::cout << std::endl;
    }
}

/*!
 *! \brief Solve a system non-linear eq. with the modified Newton
 */
void mod_newt_sys() {
    
    // Function parameters
    MatrixXd A = MatrixXd::Random(4,4);
    VectorXd c = VectorXd::Random(4);
    
    A = (A*A.transpose()).eval(); // make sure A is symmetric
    c = c.cwiseAbs(); // make sure c is > 0 uniformly
    
    // Handler for function and jacobian in standard format. Must see matrix A and vector c
    auto F = [&A,&c] (const VectorXd & x) { 
        VectorXd tmp =  A*x + c.cwiseProduct(x.array().exp().matrix()).eval();
        return tmp;
    };
    auto dF = [&A, &c] (const VectorXd & x) { 
        MatrixXd C = A; 
        VectorXd temp = c.cwiseProduct(x.array().exp().matrix()); 
        C += temp.asDiagonal(); 
        return C; 
    };
    
    // Container for errors
    std::vector<double> err;
    // Define lambda for breaking condition, which also stores 
    // the error of the previous step. Must see err
    auto rerr = [&err] (VectorXd & x) { 
        if(err.size() > 0) 
            return err.back();
        else return 1.;
        
    };
    
    std::cout << std::endl << "*** Modified Newton method (system) ***" 
              << std::endl << std::endl;
    VectorXd x_system = VectorXd::Zero(c.size()); // Initial guess
    
    // Refactor (i.e. wrap) the step st. is compatible with general\_nonlinear\_solver,
    // must see F; dF and error vector
    // We also store the error
    auto newt_system_step = [&F, &dF, &err] (const VectorXd & x) -> VectorXd {
        VectorXd x_new = mod_newt_step_system(x, F, dF);
        double e = (x - x_new).norm()/ x_new.norm();
        err.push_back(e);
        return x_new;
    };
    
    // Actually performs computations
    general_nonlinear_solver(newt_system_step, x_system, rerr);
    
    // Sanity check
    std::cout << std::endl << "x^*_system = " 
              << std::endl << x_system 
              << std::endl;
    std::cout << std::endl << "F(x^*_system) = " 
              << std::endl << F(x_system) 
              << std::endl;
}

int main() {
    // Part 1: solve scalar
    mod_newt_ord();
    
    // Part 2: solve system
    mod_newt_sys();
}
