#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/KroneckerProduct>
#include <vector>

#include "dampnewton.hpp"

/**
 * \brief Implements a Runge-Kutta implicit solver for a
 * given Butcher tableau for autonomous ODEs.
 */
/* SAM_LISTING_BEGIN_1 */
class implicitRKIntegrator {
 public:
  /**
   * \brief Constructor for the implicit RK method.
   * Performs size checks and copies A and b into internal storage.
   * \param A Matrix containing coefficents of Butcher tableau.
   * \param b Vector containing coefficients of
   * lower part of Butcher tableau.
   */
  implicitRKIntegrator(const Eigen::MatrixXd &A, const Eigen::VectorXd &b)
      : A(A), b(b), s(b.size()) {
    assert(A.cols() == A.rows() && "Matrix must be square.");
    assert(A.cols() == b.size() && "Incompatible matrix/vector size.");
  }

  /**
   * \brief Perform the solution of the ODE.
   * Solve an autonomous ODE y' = f(y), y(0) = y0, using an
   * implicit RK scheme given in the Butcher tableau provided in the
   * constructor. Performs N equidistant steps upto time T
   * with initial data y0.
   * \tparam Function type for function implementing the rhs function.
   * Must have Eigen::VectorXd operator()(Eigen::VectorXd x)
   * \tparam Function2 type for function implementing the Jacobian of f.
   * Must have Eigen::MatrixXd operator()(Eigen::VectorXd x)
   * \param f function handle for rhs in y' = f(y), e.g.
   * implemented using lambda funciton
   * \param Jf function handle for Jf, e.g.
   * implemented using lambda funciton
   * \param T final time T
   * \param y0 initial data y(0) = y0 for y' = f(y)
   * \param N number of steps to perform.
   * Step size is h = T / N. Steps are equidistant.
   * \return vector containing all steps y^n (for each n)
   * including initial and final value
   */
  template <class Function, class Jacobian>
  std::vector<Eigen::VectorXd> solve(Function &&f, Jacobian &&Jf, double T,
                                     const Eigen::VectorXd &y0,
                                     unsigned int N) const {
    // Init step size
    double h = T / N;

    // Will contain all steps, reserve memory for efficiency
    std::vector<Eigen::VectorXd> res;
    res.reserve(N + 1);

    // Store initial data
    res.push_back(y0);
    // Initialize some memory to store temporary values
    Eigen::VectorXd ytemp1 = y0;
    Eigen::VectorXd ytemp2 = y0;
    // Pointers to swap previous value
    Eigen::VectorXd *yold = &ytemp1;
    Eigen::VectorXd *ynew = &ytemp2;

    // Loop over all fixed steps
    for (unsigned int k = 0; k < N; ++k) {
      // Compute, save and swap next step
      step(f, Jf, h, *yold, *ynew);
      res.push_back(*ynew);
      std::swap(yold, ynew);
    }
    return res;
  }

 private:
  /**
   * \brief Perform a single step of the RK method for the
   * solution of the autonomous ODE
   * Compute a single explicit RK step y^{n+1} = y_n + \sum ...
   * starting from value y0 and storing next value in y1
   * \tparam Function type for function implementing the rhs.
   * Must have Eigen::VectorXd operator()(Eigen::VectorXd x)
   * \tparam Jacobian type for function implementing the Jacobian of f.
   * Must have Eigen::MatrixXd operator()(Eigen::VectorXd x)
   * \param f function handle for ths f, s.t. y' = f(y)
   * \param Jf function handle for Jf, e.g. implemented using lambda
   *function
   * \param h step size ! \param[in] y0 initial Eigen::VectorXd
   * \param y1 next step y^{n+1} = y^n + ...
   */
  template <class Function, class Jacobian>
  void step(Function &&f, Jacobian &&Jf, double h, const Eigen::VectorXd &y0,
            Eigen::VectorXd &y1) const {
    // TODO: (12-1.c) Perform one step of the implicit Runge-Kutta method.
    // START
    const unsigned int d = y0.size();
    Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(d, d);

    // Handle for the function F describing the
    // equation satisfied by the stages g
    auto F = [&y0, h, d, this, &f, &eye](Eigen::VectorXd gv) {
      Eigen::VectorXd Fv = gv;
      for (unsigned int j = 0; j < s; j++) {
        Fv = Fv -
             h * kroneckerProduct(A.col(j), eye) * f(y0 + gv.segment(j * d, d));
      }
      return Fv;
    };

    // Handle for the Jacobian of F.
    auto JF = [&y0, h, d, &Jf, this](Eigen::VectorXd gv) {
      Eigen::MatrixXd DF(s * d, s * d);
      for (unsigned int j = 0; j < s; j++) {
        DF.block(j * d, 0, d, s * d) =
            kroneckerProduct(A.row(j), Jf(y0 + gv.segment(j * d, d)));
      }
      DF = Eigen::MatrixXd::Identity(s * d, s * d) - h * DF;
      return DF;
    };

    // Obtain stages with damped Newton method
    Eigen::VectorXd gv = Eigen::VectorXd::Zero(s * d);
    dampnewton(F, JF, gv);

    // Calculate y1
    Eigen::MatrixXd K(d, s);
    for (unsigned int j = 0; j < s; j++)
      K.col(j) = f(y0 + gv.segment(j * d, d));
    y1 = y0 + h * K * b;
    // END
  }
  //<! Matrix A in Butcher scheme
  const Eigen::MatrixXd A;
  //<! Vector b in Butcher scheme
  const Eigen::VectorXd b;
  //<! Size of Butcher matrix and vector A and b
  unsigned int s;
};
/* SAM_LISTING_END_1 */
