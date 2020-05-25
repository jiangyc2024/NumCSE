#ifndef IMPLICIT_RKINTEGRATOR_HPP
#define IMPLICIT_RKINTEGRATOR_HPP

//#pragma once

#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "dampnewton.hpp"

using namespace Eigen;

/*!
 *! \file implicit_rkintegrator.hpp A header only file implementing a
 *! implicit runge kutta integrator.
 */

#ifndef KRON
#define KRON
/*!
 *! \brief Compute the Kronecker product $C = A \otimes B$.
 *! \param[in] A Matrix $m \times n$.
 *! \param[in] B Matrix $l \times k$.
 *! \return Kronecker product of A and B of dim $ml \times nk$.
 */
MatrixXd kron(const MatrixXd &A, const MatrixXd &B) {
  MatrixXd C(A.rows() * B.rows(), A.cols() * B.cols());
  for (unsigned int i = 0; i < A.rows(); ++i) {
    for (unsigned int j = 0; j < A.cols(); ++j) {
      C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
    }
  }
  return C;
}
#endif

/* SAM_LISTING_BEGIN_1 */
/*!
 *! \brief Implements a Runge-Kutta implicit solver for a
 *! given Butcher tableau for autonomous ODEs.
 */
class implicit_RKIntegrator {
public:
  /*!
   *! \brief Constructor for the implicit RK method.
   *! Performs size checks and copies A and b into internal storage.
   *! \param[in] A Matrix containing coefficents of Butcher tableau.
   *! \param[in] b Vector containing coefficients of
   *! lower part of Butcher tableau.
   */
  implicit_RKIntegrator(const MatrixXd &A, const VectorXd &b)
      : A(A), b(b), s(b.size()) {
    assert(A.cols() == A.rows() && "Matrix must be square.");
    assert(A.cols() == b.size() && "Incompatible matrix/vector size.");
  }

  /*!
   *! \brief Perform the solution of the ODE.
   *! Solve an autonomous ODE y' = f(y), y(0) = y0, using an
   *! implicit RK scheme given in the Butcher tableau provided in the
   *! constructor. Performs N equidistant steps upto time T
   *! with initial data y0.
   *! \tparam Function type for function implementing the rhs function.
   *! Must have VectorXd operator()(VectorXd x)
   *! \tparam Function2 type for function implementing the Jacobian of f.
   *! Must have MatrixXd operator()(VectorXd x)
   *! \param[in] f function handle for rhs in y' = f(y), e.g.
   *! implemented using lambda funciton
   *! \param[in] Jf function handle for Jf, e.g.
   *! implemented using lambda funciton
   *! \param[in] T final time T
   *! \param[in] y0 initial data y(0) = y0 for y' = f(y)
   *! \param[in] N number of steps to perform.
   *! Step size is h = T / N. Steps are equidistant.
   *! \return vector containing all steps y^n (for each n)
   *! including initial and final value
   */
  template <class Function, class Jacobian>
  std::vector<VectorXd> solve(Function &&f, Jacobian &&Jf, double T,
                              const VectorXd &y0, unsigned int N) const {
    // Iniz step size
    double h = T / N;

    // Will contain all steps, reserve memory for efficiency
    std::vector<VectorXd> res;
    res.reserve(N + 1);

    // Store initial data
    res.push_back(y0);
    // START
    // Your solution here
    // END
    return res;
  }

private:
  // START
  // Any auxiliary functions or data here
  // END
  //<! Matrix A in Butcher scheme
  const MatrixXd A;
  //<! Vector b in Butcher scheme
  const VectorXd b;
  //<! Size of Butcher matrix and vector A and b
  unsigned int s;
};
/* SAM_LISTING_END_1 */

#endif
