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

//#include "matplotlibcpp.h"

using namespace Eigen;
//namespace plt = matplotlibcpp;

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
 * @param[in] t $n$ size vector
 * @param[out] gamma $4$ size vector
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd data_fit_normal(const VectorXd &t, const VectorXd &b) {
  //TO DO (4-3.a) Solve normal equations to find the coefficients of the 
  // linear fitting
  VectorXd gamma(4);
  //START
  MatrixXd A = make_A(t);
  MatrixXd At = A.transpose();
  MatrixXd AtA = At * A;
  VectorXd Atb = At * b;
  gamma =  AtA.ldlt().solve(Atb);
  //END
  return gamma;
}
/* SAM_LISTING_END_2 */



/* @param[in] b $n$ size vector 
 * @param[in] t $n $ size vector
 * @param[out] gamma $4$ size vector
 */
/* SAM_LISTING_BEGIN_3 */
VectorXd data_fit_qr(const VectorXd &t, const VectorXd &b) {
  //TO DO (4-3.b) Find the coefficients for the linear
  // fitting by means of the QR decomposition of A
  VectorXd gamma(4);
  //START
  MatrixXd A = make_A(t);
  gamma = A.colPivHouseholderQr().solve(b);
  //END
  return gamma;
}
/* SAM_LISTING_END_3 */



/* @param[in] gamma $4$ size column vector
 * @param[in] t  vector 
 * @param[out] y  vector of size = t.size() 
 */
//Note: the code will not run until this function is implemented
 /* SAM_LISTING_BEGIN_4 */
VectorXd fit_plot(const VectorXd &gamma, const VectorXd &t) {
  //TO DO (4-3.c): Define the data for the first plot:
  // evaluate the function f at the grid defined by t

  VectorXd y;
  //START
  MatrixXd A = make_A(t);
  y = A*gamma;
  //END
  return y;
}
/* SAM_LISTING_END_4 */


/* @param[in] gamma $4$ size column vector 
 * @param[out] err $n$ size vector
 */
 //Note: the code will not run until this function is implemented
 /* SAM_LISTING_BEGIN_5 */
VectorXd error_plot( const VectorXd &gamma ) {
  //TO DO (4-3.c): Compute the vector of squared errors of your 
  // fit at the provided sample points
  VectorXd err;
  //START
  VectorXd t = VectorXd::LinSpaced(10, 0.1, 1.0);
  VectorXd f(10);
  f << 100. , 34. , 17. , 12. , 9. , 6. , 5. , 4. , 4. , 2.;

  err = fit_plot(gamma,t)-f;
  err = err.cwiseProduct(err);
  //END
  return err;
}
/* SAM_LISTING_END_5 */

#endif
