#ifndef RKINTEGRATOR_HPP
#define RKINTEGRATOR_HPP

#include <Eigen/Dense>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * @brief A Runge-Kutta explicit solver for a given
 * Butcher tableau for autonomous ODEs.
 *
 * @tparam State State a type representing the space in which the solution
 * lies, e.g. R^d, represented by e.g. VectorXd.
 */
/* SAM_LISTING_BEGIN_0 */
template <class State>
class RKIntegrator {
 public:
  /**
   * @brief Constructor for the RK method.
   * Performs size checks and copies $\VA$ and $\Vb$ into internal storage.
   *
   * @param A Matrix containing coefficents of the Butcher tableau,
   * must be (strictly) lower triangular (no check is done).
   * @param b Vector containing coefficients of lower
   * part of Butcher tableau.
   */
  RKIntegrator(const Eigen::MatrixXd &A, const Eigen::VectorXd &b)
  // TODO: (11-2.a) Use an initializer list to Initialize the private
  // variables A, b, s.
  // START

  // END
  {
    // TODO: (11-2.a) Use assert() to make sure A is a square matrix of
    // shape s*s with s=b.size().
    // START

    // END
  }

  /**
   * @brief Perform the solution of the ODE.
   * Solve an autonomous ODE $y' = f(y)$, $y(0) = y0$, using a
   * RK scheme given in the Butcher tableau provided in the
   * constructor. Performs $N$ equidistant steps up to time
   * $T$ with initial data $y_0$.
   *
   * @tparam Function type for function implementing the rhs function.
   * Must have State operator()(const State & x).
   * @param f function handle for rhs in $y' = f(y)$, e.g.
   * implemented using lambda funciton.
   * @param T The final time $T$.
   * @param y0 Initial data $y(0) = y_0$ for $y' = f(y)$.
   * @param N  Number of steps to perform. Step size is $h = T / N$.
   * Steps are equidistant.
   * @return std::vector<State> containing all steps $y^n$ (for each $n$)
   * including initial and final value.
   */
  template <class Function>
  std::vector<State> solve(Function &&f, double T, const State &y0,
                           unsigned int N) const {
    std::vector<State> res;
    // Will contain all steps, reserve memory for efficiency
    res.reserve(N + 1);
    // TODO: (11-2.a) Compute the solution of y'=f(y), y(0)=y0 up to time T
    // using N Runge-Kutta steps, and store the solution at step n in res[n].
    // Use the private method step() for the Runge-Kutta step.
    // START

    // END
    return res;
  }

 private:
  /**
   * @brief Perform a single step of the RK method.
   * Solve an autonomous ODE using an explicit Runge Kutta Method.
   * Compute a single explicit RK step $y^{n+1} = y_n + \sum \dots$
   * starting from value $y_0$ and storing next value in $y_1$.
   *
   * @tparam Function  type for function implementing the rhs.
   * Must have State operator()(State x)
   * @param f function handle for rhs $f$, s.t. $y' = f(y)$
   * @param h step size
   * @param y0 initial state
   * @param y1 next step $y^{n+1} = y^n + \dots$
   */
  template <class Function>
  void step(const Function &f, double h, const State &y0, State &y1) const {
    // TODO: (11-2.a) Compute a single step of the explicit Runge-Kutta method
    // (defined by A, b, s) for y'=f(y) starting at State y0 and using
    // step size h. Store the result in the State y1.
    // START

    // END
  }

  const Eigen::MatrixXd A;  //!< Matrix $\VA$ in Butcher scheme
  const Eigen::VectorXd b;  //!< Vector $\Vb$ in Butcher scheme
  unsigned int s;  //!< Size of Butcher matrix and vector $\VA$ and $\Vb$
};
/* SAM_LISTING_END_0 */

#endif