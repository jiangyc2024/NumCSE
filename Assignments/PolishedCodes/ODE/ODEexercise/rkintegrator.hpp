#ifndef RKINTEGRATOR_HPP
#define RKINTEGRATOR_HPP

#include <Eigen/Dense>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * \brief A Runge-Kutta explicit solver for a given
 * Butcher tableau for autonomous ODEs.
 *
 * \tparam State State a type representing the space in which the solution
 * lies, e.g. R^d, represented by e.g. VectorXd.
 * \tparam RHSFunctor a functor that takes State as input and return a State
 *
 */
template <class State, class RHSFunctor>
class RKIntegrator {
 public:
  /**
   * \brief Constructor for the RK method.
   * Performs size checks and copies $\VA$ and $\Vb$ into internal storage.
   *
   * \param f Right hand side function f of the ODE $\dot(y) = f(y)$
   * \param A Matrix containing coefficents of the Butcher tableau,
   * must be (strictly) lower triangular (no check is done).
   * \param b Vector containing coefficients of lower
   * part of Butcher tableau.
   */
  RKIntegrator(RHSFunctor& f, const Eigen::MatrixXd &A, const Eigen::VectorXd &b)
      : f_(f), A_(A), b_(b) {};

  /**
   * \brief Perform a single step of the RK method.
   * Solve an autonomous ODE using an explicit Runge Kutta Method.
   * Compute a single explicit RK step $y^{n+1} = y_n + \sum \dots$
   * starting from value $y_0$ and storing next value in $y_1$.
   *
   * \tparam Function  type for function implementing the rhs.
   * Must have State operator()(State x)
   * \param f function handle for rhs $f$, s.t. $y' = f(y)$
   * \param h step size
   * \param y0 initial state
   * \param y1 next step $y^{n+1} = y^n + \dots$
   */
  /* SAM_LISTING_BEGIN_1 */
  State step(double h, const State &y0) const {
    State y1{y0};

    // TODO: (0-2.e) Compute a single step of the explicit Runge-Kutta method
    // (defined by A, b, s) for y'=f(y) starting at State y0 and using
    // step size h. Store the result in the State y1.
    // START

    // Reserve space for increments
    std::vector<State> k;
    k.reserve(b_.size());

    // Loop over the size of RK
    for (unsigned int i = 0; i < b_.size(); ++i) {
      // Compute increments and save them to $k$
      State incr = y0;
      for (unsigned int j = 0; j < i; ++j) {
        incr += h * A_(i, j) * k.at(j);
      }
      k.emplace_back(f_(incr));
      y1 += h * b_(i) * k.back();
    }
    // END
    return y1;
  }
  /* SAM_LISTING_END_1 */

  /**
   * \brief Perform the solution of the ODE.
   * Solve an autonomous ODE $y' = f(y)$, $y(0) = y0$, using a
   * RK scheme given in the Butcher tableau provided in the
   * constructor. Performs $N$ equidistant steps up to time
   * $T$ with initial data $y_0$.
   *
   * \tparam Function type for function implementing the rhs function.
   * Must have State operator()(const State & x).
   * \param f function handle for rhs in $y' = f(y)$, e.g.
   * implemented using lambda funciton.
   * \param T The final time $T$.
   * \param y0 Initial data $y(0) = y_0$ for $y' = f(y)$.
   * \param N  Number of steps to perform. Step size is $h = T / N$.
   * Steps are equidistant.
   * \return std::vector<State> containing all steps $y^n$ (for each $n$)
   * including initial and final value.
   */
   /* SAM_LISTING_BEGIN_2 */
  std::vector<State> solve(double T, const State &y0, unsigned int N) const {
    std::vector<State> ys; // contain states at all steps

    // TODO: (0-2.f) Compute the solution of y'=f(y), y(0)=y0 up to time T
    // using N Runge-Kutta steps, and store the solution at step n in ys[n].
    // Use the member method step() for the Runge-Kutta step.

    // START
    // Step size
    const double h = T / N;
    // reserve memory for efficiency
    ys.reserve(N + 1);
    // Store initial data
    ys.push_back(y0);

    State y = y0;
    // Loop over all fixed steps
    for (unsigned int k = 0; k < N; ++k) {
      // Compute, save and swap next step
      y = step(h, y);
      ys.push_back(y);
    }
    // END
    return ys;
  }
  /* SAM_LISTING_END_2 */

 private:

  RHSFunctor f_;            //!< rhs function $f$
  const Eigen::MatrixXd A_;  //!< Matrix $\VA$ in Butcher scheme
  const Eigen::VectorXd b_;  //!< Vector $\Vb$ in Butcher scheme
};

/* SAM_LISTING_BEGIN_3 */
void testcvgRK() {
  // Matrices involved in the ODE system
  Eigen::Matrix2d M1;
  M1 << 1., 2., 1., 0;
  Eigen::Matrix2d M2;
  M2 << -1., 0, 2., 2.;

  // Implementation of butcher scheme (0.2.1)
  Eigen::Matrix2d A;
  Eigen::Vector2d b;
  A << 0, 0, 2./3., 0;
  b << 1./4., 3./4.;

  // TODO: (0-2.h) Test the convergence rate of RK scheme (0.2.1) when solving
  // the second-order ODE (0.2.2). You can directly invoke solve() function assuming that it works.
  // Hint:
  //    Step1. compute a reference solution with sufficiently large N (e.g. 4096).
  //    Step2. compute a series of solution under increasing N (e.g. 16, 32, 64, ...,1024). Output the error.
  //    Step3. conclude the convergence rate by inspecting the errors under different N.
  //
  // START

  // Notice that the 1st-order ODE we are actually solving is of dimension 4.
  // Construct right hand side function
  auto f = [&M1, &M2](const Eigen::Vector4d& y)->Eigen::Vector4d {
    Eigen::Matrix4d mat;
    mat.block(0,0,2,2) = Eigen::Matrix2d::Zero();
    mat.block(0,2,2,2) = Eigen::Matrix2d::Identity();
    mat.block(2,0,2,2) = - M2;
    mat.block(2,2,2,2) = - M1;
    return mat * y;
  };

  // Initialize RK with Butcher scheme and the rhs function
  RKIntegrator<Eigen::Vector4d, decltype(f)> RK(f, A, b);

  // Initial value for model. Be careful about the order.
  Eigen::Vector4d y0;
  y0 << 1., 0, 2., 3.;

  // Final time for model
  constexpr double T = 1.;

  auto ref_sol = RK.solve(T, y0, 4096).back();

  // Array of number of steps
  Eigen::ArrayXd N(6);
  N << 32, 64, 128, 256, 512, 1024;
  Eigen::ArrayXd Error(N.size());

  // Start convergence study
  std::cout << std::setw(15) << "N" << std::setw(15) << "error" << std::endl;
  for (unsigned int i = 0; i < N.size(); ++i) {
    auto sol = RK.solve(T, y0, N[i]).back();
    const double err = (sol - ref_sol).norm();
    std::cout << std::setw(15) << N[i] << std::setw(15) << err << std::endl;
  }

  // END
  return;
}
/* SAM_LISTING_END_3 */

#endif
