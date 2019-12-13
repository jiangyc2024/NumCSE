#ifndef ODESOLVE_HPP
#define ODESOLVE_HPP 

#include <iostream>
#include <vector>


#include <Eigen/Dense>

#include "matplotlibcpp.h"
#include "polyfit.hpp"

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

namespace plt = matplotlibcpp;


//! \tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
//! \param[in] Psi original evolution operator, must have operator(double, const Vector&)
//! \param[in] p parameter p for construction of Psi tilde
//! \param[in] h step-size
//! \param[in] y0 previous step
/* SAM_LISTING_BEGIN_0 */
template <class DiscEvlOp>
Vector psitilde(DiscEvlOp &&Psi, unsigned int p, double h, const Vector &y0) {
  // TO DO (12-3.b): apply the evolution operator \tilde{\Psi}
  // with step-size h to the value y0
  // START
  
  // END
  return y0;
}
/* SAM_LISTING_END_0 */


//! \tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
//! \param[in] Psi original evolution operator, must have operator(double, const Vector&)
//! \param[in] T final time
//! \param[in] y0 initial data
//! \param[in] N number of steps
/* SAM_LISTING_BEGIN_1 */
template <class DiscEvlOp>
std::vector<Vector> odeintequi(DiscEvlOp &&Psi, double T, const Vector &y0, unsigned int N) {
  std::vector<Vector> Y;
  // TO DO (12-3.c): Compute y from time 0 to T using N equidistant time steps
  // return a std::vector containing all steps y_0,...,y_N
  // START
  
  // END
  return Y;
}
/* SAM_LISTING_END_1 */



/* SAM_LISTING_BEGIN_2 */
double testcvpExtrapolatedEuler(void) {
  
  double conv_rate;
  double T = 1.;
  Vector y0 = Vector::Zero(1);
  auto f = [] (const Vector &y) -> Vector { return Vector::Ones(1) + y*y; };
  // TO DO (12-3.d): tabulate the values of the error corresponding to 
  // \tilde{\psi}, where \psi is the explicit Euler method. 
  // return the empirical convergence rate using polyfit.
  // Hint: first define a lambda for \psi. Then use psitilde to obtain a
  // suitable input for odeintequi.
  
  // START
  
  // END
  return conv_rate;
}
/* SAM_LISTING_END_2 */




//! \tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
//! \param[in] Psi low-order evolution operator, must have operator(double, const Vector&)
//! \param[in] T final time
//! \param[in] y0 initial data
//! \param[in] h0 initial step size
//! \param[in] p parameter p for construction of Psi tilde
//! \param[in] reltol relative tolerance for error control
//! \param[in] abstol absolute tolerance for error control
//! \param[in] hmin minimal step size
/* SAM_LISTING_BEGIN_3 */
template <class DiscEvlOp>
std::pair< std::vector<double>, std::vector<Vector> > 
odeintssctrl(DiscEvlOp&& Psi, double T, const Vector &y0, double h0, 
             unsigned int p, double reltol, double abstol, double hmin) {
  std::vector<double> t;
  std::vector<Vector> Y;               
  // TO DO (12-3.e):  Compute y from time 0 to T with adaptive time step.
  // Display a warning if the tolerance cannot be met with minimum  
  // step size hmin. return a pair of vectors containing the times and
  // the computed values.
  // START
  
  // END
  return std::make_pair(t,Y);
}
/* SAM_LISTING_END_3 */


/* SAM_LISTING_BEGIN_4 */
void solveTangentIVP(void)
{
  auto f = [] (const Vector &y) -> Vector { return Vector::Ones(1) + y*y; };
  Vector y0 = Vector::Zero(1);
  // TO DO (12-3.f): run the adaptive integration algorithm and plot the 
  // resulting values of y(t). 
  // Hint: you might use a loop to convert a std::vector<Vector> into a
  // std::vector<double>, since each Vector has size 1 
  
  // START
  
  // END
}
/* SAM_LISTING_END_4 */
#endif