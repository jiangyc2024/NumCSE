#pragma once

#include <vector>
#include <cassert>

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

//! \file rkintegrator.hpp Implementation of RkIntegrator class.

/*!
 *! @brief A Runge-Kutta explicit solver for a given
 *! Butcher tableau for autonomous ODEs.
 *! @tparam State a type representing the space in which the solution
 *! lies, e.g. R^d, represented by e.g. Eigen::VectorXd.
 */
/* SAM_LISTING_BEGIN_0 */
template <class State>
class RKIntegrator {
public:
    /*!
     *! @brief Constructor for the RK method.
     *! Performs size checks and copies A and b into internal storage.
     *! @param[in] A Matrix containing coefficents of the Butcher tableau,
     *! must be (strictly) lower triangular (no check is done).
     *! @param[in] b Vector containing coefficients of lower
     *! part of Butcher tableau.
     */
    RKIntegrator(const Eigen::MatrixXd & A, const Eigen::VectorXd & b)
        : A(A), b(b), s(b.size()) {
        assert( A.cols() == A.rows() && "Matrix must be square.");
        assert( A.cols() == b.size() && "Incompatible matrix/vector size.");
    }

    /*!
     *! @brief Perform the solution of the ODE
     *! Solve an autonomous ODE y' = f(y), y(0) = y0, using a RK scheme given in the Butcher tableau provided in the
     *! constructor. Performs N equidistant steps upto time T with initial data y0
     *! @tparam Function type for function implementing the rhs function. Must have State operator()(State x)
     *! @param[in] f function handle for rhs in y' = f(y), e.g. implemented using lambda funciton
     *! @param[in] T final time T
     *! @param[in] y0 initial data y(0) = y0 for y' = f(y)
     *! @param[in] N number of steps to perform. Step size is h = T / N. Steps are equidistant.
     *! @return vector containing all steps y^n (for each n) including initial and final value
     */
    template <class Function>
    std::vector<State> solve(const Function &f, double T, const State & y0, unsigned int N) const {
        std::vector<State> res;
        // TODO: computes $N$ uniform time steps for the ODE $y'(t) =
        // f(y)$ up to time $T$ of RK method with initial value $y0$ and store
        // all steps $y_k$ into return vector
        return res;
    }
private:
    /*!
     *! @brief Perform a single step of the RK method.
     *! Solve an autonomous ODE using an explicit Runge Kutta Method.
     *! Compute a single explicit RK step y^{n+1} = y_n + \sum ...
     *! starting from value y0 and storing next value in y1.
     *! @tparam Function type for function implementing the rhs. Must have State operator()(State x)
     *! @param[in] f function handle for ths f, s.t. y' = f(y)
     *! @param[in] h step size
     *! @param[in] y0 initial state
     *! @param[out] y1 next step y^{n+1} = y^n + ...
     */
    template <class Function>
    void step(const Function &f, double h,
              const State & y0, State & y1) const {
        // TODO: performs a single step from $y0$ to $y1$ with
        // step size $h$ of the RK method for the IVP with rhs $f$
    }

    //!< Matrix A in Butcher scheme
    const Eigen::MatrixXd A;
    //!< Vector b in Butcher scheme
    const Eigen::VectorXd b;
    //!< Size of Butcher matrix and vector A and b
    unsigned int s;
};
/* SAM_LISTING_END_0 */
