/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */
#ifndef PERIODICCOLLOCATIONHPP
#define PERIODICCOLLOCATIONHPP

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <unsupported/Eigen/FFT>
#include <vector>
#include "helper_functions.hpp"

// Evaluation of a trigonometric polynomials in M equidistant points $\cob{M\geq
// N}$
Eigen::VectorXd eval_uN(const Eigen::VectorXd &x, unsigned int M) {
  unsigned int N = x.size() - 1;
  assert((M > N) && "Number of evaluation points must be > N");
  Eigen::VectorXd u(M);
  // TO DO: Write efficient implementation to evaluate u_N using coefficients in
  //        x at M equidistant points 
  // START
  
  // END
  return u;
}

Eigen::VectorXd eval_F(const Eigen::VectorXd &x) {
  const unsigned int N = x.size() - 1;
  Eigen::VectorXd Fx(N+1);
  // TO DO: Write a function to evaluate F using eval_uN
  // START 
 
  // END
  return Fx;
}

Eigen::MatrixXd eval_DF(const Eigen::VectorXd &x) {
  const unsigned int N = x.size() - 1;
  Eigen::MatrixXd J(N + 1, N + 1);
  // TO DO: Write a function to evaluate the Jacobian of F
  // START
  
  //END
  return J;
}

#endif
