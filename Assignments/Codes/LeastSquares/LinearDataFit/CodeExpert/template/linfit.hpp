#ifndef LINFIT_HPP
#define LINFIT_HPP

//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <cmath>

#include <Eigen/Dense>

#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;

/* @param[in] b $n$ size vector 
 * @param[out] A $n \times 4$ matrix
 */
/* SAM_LISTING_BEGIN_1 */
MatrixXd make_A(const VectorXd &b) {
  size_t n = b.size();
  MatrixXd A(n, 4);
  //TO DO (4-3.a) Build the matrix A 
  // Hint: evaluate the functions \phi_j at the time points defined in b
  //START
  for (size_t i = 0; i < n; i++) {
    double t = b(i);
    A(i, 0) = 1.0 / t;
    A(i, 1) = 1.0 / (t*t);
    A(i, 2) = std::exp(-(t-1));
    A(i, 3) = std::exp(-2*(t-1));
  }
  //END 
  return A;
}
/* SAM_LISTING_END_1 */



/* @param[in] b $n$ size vector 
 * @param[out] A $n \times 4$ matrix
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd data_fit_normal(const MatrixXd &A, const VectorXd &b) {
  //TO DO (4-3.a) Solve normal equations to find the coefficients of the 
  // linear fitting
  //START
  MatrixXd At = A.transpose();
  MatrixXd AtA = At * A;
  VectorXd Atb = At * b;
  return AtA.ldlt().solve(Atb);
  //END
}
/* SAM_LISTING_END_2 */



/* @param[in] b $n$ size vector 
 * @param[out] A $n \times 4$ matrix
 */
/* SAM_LISTING_BEGIN_3 */
VectorXd data_fit_qr(const MatrixXd &A, const VectorXd &b) {
  //TO DO (4-3.b) Find the coefficients for the linear
  // fitting by means of the QR decomposition of A
  //START
  return A.colPivHouseholderQr().solve(b);
  //END
}
/* SAM_LISTING_END_3 */



/* @param[in] gamma $4$ size column vector
 * @param[in] tl high-resolution vector 
 * @param[out] yl1 high-resolution vector 
 * @param[out] yl2 high-resolution vector 
 */
 /* SAM_LISTING_BEGIN_4 */
void fit_plot(const VectorXd &gamma, const VectorXd &tl,
                VectorXd &yl) {
  //TO DO (4-3.c): Define the data for the first plot:
  // evaluate the function f at the high-resolution grid defined by tl
  //START
  MatrixXd Al = make_A(tl);
  yl =  Al*gamma;
  //END
}
/* SAM_LISTING_END_4 */


/* @param[in] gamma1 $4$ size column vector 
 * @param[in] A $n \times 4$ size matrix 
 * @param[in] f $n$ size vector
 * @param[out] err1 $n$ size vector
 * @param[out] err2 $n$ size vector
 */
 /* SAM_LISTING_BEGIN_5 */
void error_plot( const VectorXd &gamma, const MatrixXd &A, 
                  const VectorXd &f, VectorXd &err ) {
  //TO DO (4-3.c): Define the data for the second plot:
  // evaluate the function at the data points and compute the
  // l^2 error squared. Here, the vector f contains the measured samples. 
  //START
  VectorXd y = A*gamma;
  // compute squared errors at each data point
  err = y-f;
  err = err.cwiseProduct(err);
  //END
}
/* SAM_LISTING_END_5 */

#endif