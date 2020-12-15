#ifndef MODNEWT_HPP
#define MODNEWT_HPP

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

double norm(double x) { return std::abs(x); }
double norm(const Eigen::VectorXd& x) { return x.norm(); }

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
/* SAM_LISTING_BEGIN_1 */
template <typename Scalar, class Function, class Jacobian>
Scalar mod_newt_step_scalar(const Scalar& x, Function&& f, Jacobian&& df) {
  double return_value = 0.;
  // TODO: (9-9.b) Compute a step of the modified Newton method for a scalar
  // function START

  // END
  return return_value;
}
/* SAM_LISTING_END_1 */

/*!
 *! \brief Performs many steps of an iteration and terminate when convergence
 *reached ! or maximum number of iterations has been reached. ! \tparam
 *StepFunction type for the step function handle ! \tparam Vector argument type
 *passed to the iteration function ! \tparam ErrorFunction type for the error
 *function computing error of the method ! \param[in] step Function implementing
 *the step x_{k+1} = step(x_{k}), signature Vector(const Vector&) !
 *\param[in,out] x initial data (as input) and final iteration (as output) !
 *\param[in] errf function implementing the norm of the error (errf(x)) for
 *termination condition ! \param[in] eps tolerance to break iterations when
 *res(x) < eps ! \param[in] max_itr maximal number of iterations
 */
/* SAM_LISTING_BEGIN_3 */
template <class StepFunction, class Vector, class ErrorFunction>
bool sample_nonlinear_solver(const StepFunction& step, Vector& x,
                             const ErrorFunction& errf, double eps = 1e-8,
                             int max_itr = 100) {
  // TODO: (9-9.c,f) Implement this generic function that iteratively solves a
  // nonlinear equation. The function should return whether the iteration
  // converged in the maximum number of iterations. START

  // END
  return false;
}
/* SAM_LISTING_END_3 */

/**
 *! \brief Solve a scalar non-linear eq. with the modified Newton.
 * Test the empirical order of convergence of the method.
 */
/* SAM_LISTING_BEGIN_2 */
void mod_newt_ord(void) {
  // Setting up values, functions and jacobian
  constexpr double a = 0.123;
  auto f = [](double x) {
    return atan(x) - a;
  };                    // Function which sees a as it is constexpr
  double x_scalar = 5;  // initial guess

  // TODO: (9-9.c) Generate suitable terminal output to determine the order of
  // convergence of the modified Newton method applied to the scalar equation
  // that is given by $f = 0$. You may want to implement the templated function
  // sample\_nonlinear\_solver for the iterations. START

  // END
}
/* SAM_LISTING_END_2 */

/*!
 *! \brief Implements a single step of the modified newton
 *! \tparam Vector type of argument to function f: such as double or vector
 *etc... ! \tparam Function type for the function f, likely a lambda function !
 *\tparam Jacobian type for the jacobian df, likely a lambda function !
 *\param[in] x previous value to use in step, also initial guess if needed !
 *\param[in] f function handle for f(x) = 0 ! \param[in] df function handle for
 *jacobian df of f ! \return x_next next step x_{k+1} of modified Newton
 */
/* SAM_LISTING_BEGIN_4 */
template <typename Vector, class Function, class Jacobian>
Vector mod_newt_step_system(const Vector& x, Function&& f, Jacobian& df) {
  Vector return_value = x;  // this is just a dummy, you should overwrite it
  // TODO: (9-9.d) Efficiently perform one step of the modified Newton method
  // for a system of equations. START

  // END
  return return_value;
}
/* SAM_LISTING_END_4 */

/*!
 *! \brief Solve a system non-linear eq. with the modified Newton
 *! \param[in] A coefficient matrix
 *! \param[in] c coefficients of the exponential term
 *! \param[in] tol tolerance to break iterations
 *! \param[in] maxit maximum number of iterations
 *! \return solution of the non-linear system of equations
 */
/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd mod_newt_sys(const Eigen::MatrixXd& A, const Eigen::VectorXd& c,
                             double tol = 1.0e-6, int maxit = 100) {
  Eigen::VectorXd x_system = Eigen::VectorXd::Zero(c.size());  // Initial guess

  // TODO: (9-9.f) Solve the non-linear system of equations (9.9.6) using the
  // modified Newton method. You may want to use your implementation of
  // sample\_nonlinear\_solver. START

  // END
  return x_system;
}
/* SAM_LISTING_END_5 */

/**
 *! \brief Solve a system of non-linear eq. with the modified Newton.
 * Test the empirical order of convergence of the method.
 */
/* SAM_LISTING_BEGIN_6 */
void mod_newt_sys_test(void) {
  Eigen::MatrixXd A(4, 4);
  A << 2., -1., 0., 0., -1., 2., -1., 0., 0., -1., 2., -1., 0., 0., -1., 2.;
  Eigen::VectorXd c(4);
  c << 1., 2., 3., 4.;

  // TODO: (9-9.g) Determine the order of convergence from suitable tabulated
  // values. START

  // END
}
/* SAM_LISTING_END_6 */

#endif
