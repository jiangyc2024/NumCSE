#pragma once

#include <vector>
#include <cassert>

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

using namespace Eigen;

//! \file rkintegrator.hpp Implementation of RkIntegrator class.

/*!
 *! \brief A Runge-Kutta explicit solver for a given
 *! Butcher tableau for autonomous ODEs.
 *! \tparam State a type representing the space in which the solution
 *! lies, e.g. R^d, represented by e.g. VectorXd.
 */
/* SAM_LISTING_BEGIN_0 */
template <class State>
class RKIntegrator {
public:
  /*!
  *! \brief Constructor for the RK method.
  *! Performs size checks and copies $\VA$ and $\Vb$ into internal storage.
  *! \param[in] $\VA$ Matrix containing coefficents of the Butcher tableau,
  *! must be (strictly) lower triangular (no check is done).
  *! \param[in] $\Vb$ Vector containing coefficients of lower
  *! part of Butcher tableau.
  */
  RKIntegrator(const MatrixXd & A,
               const VectorXd & b)
    // TO DO: Use an initializer list to Initialize the private
    // variables A, b, s.
    // START
    
    // END
  {
    // TO DO: Use assert() to make sure A is a square matrix of
    // shape s*s with s=b.size().
    // START
    
    // END
  }
  
  /*!
  *! \brief Perform the solution of the ODE.
  *! Solve an autonomous ODE $y' = f(y)$, $y(0) = y0$, using a
  *! RK scheme given in the Butcher tableau provided in the
  *! constructor. Performs $N$ equidistant steps up to time
  *! $T$ with initial data $y_0$.
  *! \tparam Function type for function implementing the rhs function.
  *! Must have State operator()(const State & x).
  *! \param[in] $f$ function handle for rhs in $y' = f(y)$, e.g.\
  *! implemented using lambda funciton.
  *! \param[in] $T$ The final time $T$.
  *! \param[in] $y_0$ Initial data $y(0) = y_0$ for $y' = f(y)$.
  *! \param[in] $N$ Number of steps to perform. Step size is $h = T / N$.
  *! Steps are equidistant.
  *! \return The vector containing all steps $y^n$ (for each $n$)
  *! including initial and final value.
  */
  template <class Function>
  std::vector<State> solve(Function &&f, double T,
                           const State &y0, unsigned int N) const {
    std::vector<State> res;
    // Will contain all steps, reserve memory for efficiency
    res.reserve(N+1);
    // TO DO: Compute the solution of y'=f(y), y(0)=y0 up to time T using
    // N Runge-Kutta steps, and store the solution at step n in res[n].
    // Use the private method step() for the Runge-Kutta step.
    // START
    
    // END
    return res;
  }
  
private:
  /*!
  *! \brief Perform a single step of the RK method.
  *! Solve an autonomous ODE using an explicit Runge Kutta Method.
  *! Compute a single explicit RK step $y^{n+1} = y_n + \sum \dots$
  *! starting from value $y_0$ and storing next value in $y_1$.
  *! \tparam Function type for function implementing the rhs.
  *! Must have State operator()(State x)
  *! \param[in] $f$ function handle for rhs $f$, s.t. $y' = f(y)$
  *! \param[in] $h$ step size
  *! \param[in] $y_0$ initial state
  *! \param[out] $y_1$ next step $y^{n+1} = y^n + \dots$
  */
  template <class Function>
  void step(const Function &f, double h, const State & y0, State & y1) const {
    // TO DO: Compute a single step of the explicit Runge-Kutta method
    // (defined by A, b, s) for y'=f(y) starting at State y0 and using
    // step size h. Store the result in the State y1.
    // START
    
    // END
  }
  
  //! Matrix $\VA$ in Butcher scheme
  const MatrixXd A;
  //! Vector $\Vb$ in Butcher scheme
  const VectorXd b;
  //! Size of Butcher matrix and vector $\VA$ and $\Vb$
  unsigned int s;
};
/* SAM_LISTING_END_0 */
