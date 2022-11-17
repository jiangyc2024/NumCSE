#ifndef MODNEWT_HPP
#define MODNEWT_HPP

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

double norm(double x) { return std::abs(x); }
double norm(const Eigen::VectorXd& x) { return x.norm(); }

/**
 * @brief Implements a single step of the modified newton
 *
 * @tparam Scalar type of argument to function f: such as double, etc...
 * @tparam Function type for the function f, likely a lambda function
 * @tparam Jacobian type for the jacobian df, likely a lambda function
 * @param x previous value to use in step, also initial guess if needed
 * @param f function handle for f(x) = 0
 * @param df function handle for jacobian df of f
 * @return Scalar next step x_{k+1} of modified Newton
 */
/* SAM_LISTING_BEGIN_1 */
template <typename Scalar, class Function, class Jacobian>
Scalar mod_newt_step_scalar(const Scalar& x, Function&& f, Jacobian&& df) {
  double return_value = 0.;
  // TODO: (8-9.b) Compute a step of the modified Newton method for a scalar
  // function
  // START
  Scalar y = x + f(x) / df(x);
  return_value = y - f(y) / df(x);
  // END
  return return_value;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Performs many steps of an iteration and terminate when convergence
 * reached or maximum number of iterations has been reached.
 *
 * @tparam StepFunction type for the step function handle
 * @tparam Vector argument type passed to the iteration function
 * @tparam ErrorFunction type for the error function computing error of the
 * method
 * @param step Function implementing the step x_{k+1} = step(x_{k}),
 * signature Vector(const Vector&)
 * @param x initial data (as input) and final iteration (as output)
 * @param errf function implementing the norm of the error (errf(x)) for
 * termination condition
 * @param eps tolerance to break iterations when
 * res(x) < eps
 * @param max_itr maximal number of iterations
 * @return true if converged
 * @return false otherwise
 */
/* SAM_LISTING_BEGIN_3 */
template <class StepFunction, class Vector, class ErrorFunction>
bool sample_nonlinear_solver(const StepFunction& step, Vector& x,
                             const ErrorFunction& errf, double eps = 1e-8,
                             unsigned int max_itr = 100) {
  // TODO: (8-9.c,f) Implement this generic function that iteratively solves a
  // nonlinear equation. The function should return whether the iteration
  // converged in the maximum number of iterations
  // START.
  // Temporary where to store new step
  Vector x_new = x;
  double r = 1;

  for (unsigned itr = 0; itr < max_itr; ++itr) {
    // Compute error (or residual)
    r = errf(x);

    std::cout << "[Step " << itr << "] Error: " << r << std::endl;

    // Advance to next step, $x_{new}$ becomes $x_{k+1}$
    x_new = step(x);

    // Termination conditions
    // If tol reached, the we have convergence
    if (r < eps * norm(x)) {
      std::cout << "[CONVERGED] in " << itr << " it. due to err. err = " << r
                << " < " << eps << "." << std::endl;
      return true;
    }
    x = x_new;
  }

  // If max it reached
  std::cout << "[NOT CONVERGED] due to MAX it. = " << max_itr
            << " reached, err = " << r << "." << std::endl;
  // END
  return false;
}
/* SAM_LISTING_END_3 */

/**
 * @brief Solve a scalar non-linear eq. with the modified Newton.
 * Test the empirical order of convergence of the method.
 *
 */
/* SAM_LISTING_BEGIN_2 */
void mod_newt_ord() {
  // Setting up values, functions and jacobian
  constexpr double a = 0.123;
  auto f = [](double x) {
    return atan(x) - a;
  };                    // Function which sees a as it is constexpr
  double x_scalar = 5;  // initial guess

  // TODO: (8-9.c) Generate suitable terminal output to determine the order of
  // convergence of the modified Newton method applied to the scalar equation
  // that is given by $f = 0$. You may want to implement the templated function
  // sample\_nonlinear\_solver for the iterations.
  // START
  auto df = [](double x) { return 1. / (x * x + 1.); };  // Its derivative

  const double x_ex = tan(a);  // Exact solution

  // Store solution and error at each step
  std::vector<double> sol;
  std::vector<double> err;

  // Compute error and push back to err, used in
  // sample\_nonlinear\_solver as breaking condition errf(x) < eps
  auto errf = [&err, x_ex](double& x) {
    double e = std::abs(x - x_ex);
    err.push_back(e);
    return e;
  };

  // Perform convergence study with Modified newton for scalar
  std::cout << std::endl
            << "*** Modified Newton method (scalar) ***" << std::endl
            << std::endl;
  std::cout << "Exact: " << x_ex << std::endl;

  // Initial guess and compute initial error
  sol.push_back(x_scalar);
  errf(x_scalar);

  // Lambda performing the next step, used to define a proper
  // function handle to be passed to sample\_nonlinear\_solver
  auto newt_scalar_step = [&sol, &f, &df](double x) -> double {
    double x_new = mod_newt_step_scalar(x, f, df);
    sol.push_back(x_new);
    return x_new;
  };

  // Actually perform the solution
  sample_nonlinear_solver(newt_scalar_step, x_scalar, errf);

  // Print solution (final)
  std::cout << std::endl << "x^*_scalar = " << x_scalar << std::endl;

  // Print table of solutions, errors and EOOC
  auto space = std::setw(15);
  std::cout << space << "sol." << space << "err." << space << "order"
            << std::endl;
  for (unsigned i = 0; i < sol.size(); ++i) {
    std::cout << space << sol.at(i) << space << err.at(i);
    if (i >= 3) {
      std::cout << space
                << (log(err.at(i)) - log(err.at(i - 1))) /
                       (log(err.at(i - 1)) - log(err.at(i - 2)));
    }
    std::cout << std::endl;
  }
  // END
}
/* SAM_LISTING_END_2 */

/**
 * @brief Implements a single step of the modified newton
 *
 * @tparam Vector type of argument to function f: such as double or vector
 * etc...
 * @tparam Function type for the function f, likely a lambda function
 * @tparam Jacobian type for the jacobian df, likely a lambda function
 * @param x previous value to use in step, also initial guess if needed
 * @param f function handle for f(x) = 0
 * @param df function handle for
 * jacobian df of f
 * @return Vector next step x_{k+1} of modified Newton
 */
/* SAM_LISTING_BEGIN_4 */
template <typename Vector, class Function, class Jacobian>
Vector mod_newt_step_system(const Vector& x, Function&& f, Jacobian& df) {
  Vector return_value = x;  // this is just a dummy, you should overwrite it
  // TODO: (8-9.d) Efficiently perform one step of the modified Newton method
  // for a system of equations.
  // START
  auto lu = df(x).lu();
  // Reusing LU decomposition
  Vector y = x + lu.solve(f(x));
  return_value = y - lu.solve(f(y));
  // END
  return return_value;
}
/* SAM_LISTING_END_4 */

/**
 * @brief Solve a system non-linear eq. with the modified Newton
 *
 * @param A coefficient matrix
 * @param c coefficients of the exponential term
 * @param tol tolerance to break iterations
 * @param maxit maximum number of iterations
 * @return Eigen::VectorXd solution of the non-linear system of equations
 */
/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd mod_newt_sys(const Eigen::MatrixXd& A, const Eigen::VectorXd& c,
                             double tol = 1.0e-6, int maxit = 100) {
  Eigen::VectorXd x_system = Eigen::VectorXd::Zero(c.size());  // Initial guess

  // TODO: (8-9.f) Solve the non-linear system of equations (9.9.6) using the
  // modified Newton method. You may want to use your implementation of
  // sample\_nonlinear\_solver.
  // START

  // Handler for function and jacobian in standard format. Must see matrix $\VA$
  // and vector $\Vc$.
  auto F = [&A, &c](const Eigen::VectorXd& x) {
    Eigen::VectorXd tmp =
        A * x + c.cwiseProduct(x.array().exp().matrix()).eval();
    return tmp;
  };
  auto dF = [&A, &c](const Eigen::VectorXd& x) {
    Eigen::MatrixXd C = A;
    Eigen::VectorXd temp = c.cwiseProduct(x.array().exp().matrix());
    C += temp.asDiagonal();
    return C;
  };

  // Container for errors
  std::vector<double> err;
  // Define lambda for breaking condition, which also stores
  // the error of the previous step. Must see 'err'.
  auto rerr = [&err](Eigen::VectorXd& x) {
    if (err.size() > 0)
      return err.back();
    else
      return 1.;
  };

  // Refactor (i.e.\ wrap) the step such that it is compatible with
  // sample\_nonlinear\_solver. Must see $F$, $dF$ and 'err'. We also store the
  // error.
  auto newt_system_step = [&F, &dF,
                           &err](const Eigen::VectorXd& x) -> Eigen::VectorXd {
    Eigen::VectorXd x_new = mod_newt_step_system(x, F, dF);
    double e = (x - x_new).norm() / x_new.norm();
    err.push_back(e);
    return x_new;
  };

  // Actually performs computations
  sample_nonlinear_solver(newt_system_step, x_system, rerr, tol, maxit);
  // END
  return x_system;
}
/* SAM_LISTING_END_5 */

/**
 * @brief Solve a system of non-linear eq. with the modified Newton.
 * Test the empirical order of convergence of the method.
 *
 */
/* SAM_LISTING_BEGIN_6 */
void mod_newt_sys_test() {
  Eigen::MatrixXd A(4, 4);
  A << 2., -1., 0., 0., -1., 2., -1., 0., 0., -1., 2., -1., 0., 0., -1., 2.;
  Eigen::VectorXd c(4);
  c << 1., 2., 3., 4.;

  // TODO: (8-9.g) Determine the order of convergence from suitable tabulated
  // values.
  // determine a reference solution with very low tolerance
  // (machine precision)
  // START
  constexpr double reference_tol = std::numeric_limits<double>::epsilon();
  Eigen::VectorXd reference_sol = mod_newt_sys(A, c, reference_tol);

  std::vector<double> err;
  std::vector<Eigen::VectorXd> sol;

  // seperate calculation and printing to not pollute the table from output of
  // iteration function increase the number of maximal iterations
  sol.push_back(Eigen::VectorXd::Zero(4));
  err.push_back(reference_sol.norm());
  for (unsigned k = 1; k < 5; ++k) {
    sol.push_back(mod_newt_sys(A, c, 0., k));
    err.push_back((reference_sol - sol.at(k)).norm());
  }

  // Perform convergence study with Modified newton for system of equations
  std::cout << std::endl
            << "*** Modified Newton method (system) ***" << std::endl
            << std::endl;
  std::cout << "Exact: " << reference_sol.transpose() << std::endl;

  auto space = std::setw(15);
  std::cout << space << "k" << space << "err." << space << "order" << std::endl;

  // print iteration number, error and order
  for (unsigned k = 0; k < 5; ++k) {
    std::cout << space << k << space << err.at(k);
    if (k >= 1 && k < 4) {
      std::cout << space
                << (log(err.at(k + 1)) - log(err.at(k))) /
                       (log(err.at(k)) - log(err.at(k - 1)));
    }
    std::cout << std::endl;
  }
  // END
}
/* SAM_LISTING_END_6 */

#endif
