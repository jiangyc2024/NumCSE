#ifndef QSPLINES_HPP
#define QSPLINES_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include <vector>
#include <math.h>
#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;


/* [input] t  Vector of size n-1
   [output] pair of vectors of size n+2 and n+1, respectively
*/
/* SAM_LISTING_BEGIN_0 */
std::pair<VectorXd, VectorXd> increments (const VectorXd &t) {
  // TO DO (6-6.g) : compute the increments t_j - t_{j-1} 
  // for j = 0,...,n+1 and t_{j+1} - t_{j-1} for j = 0,...,n
  // Assume that t is sorted and that t does not include the endpoints.
  // Hint: use periodic conditions to define t_{-1}, t_{n+1} for vectorization.
  VectorXd dt;
  VectorXd ddt;
  //START
  
  // Temporary values so tests work. Replace with your own solution.
  VectorXd t1 = VectorXd::Zero(t.size() + 3);
  VectorXd t2 = VectorXd::Zero(t.size() + 2);
  return std::make_pair(t1, t2);
  
  //END
  return std::make_pair(dt,ddt);
}
/* SAM_LISTING_END_0 */




/* [input] t  Vector of size n-1
   [input] y  Vector of size n
   [output] c Vector of size n
*/
/* SAM_LISTING_BEGIN_1 */
VectorXd compute_c (const VectorXd &t, const VectorXd &y) {
  const unsigned int n = y.size(); //number of intervals
  VectorXd c(n);
  // TO DO (6-6.g) : Build the (sparse) matrix for spline interpolation 
  // at midpoints. Then compute the coefficients c for the data y. 
  
  //START
  
  //END
  return c;
}  
/* SAM_LISTING_END_1 */
  
  
  
  
  
/* [input] c  Vector of size n
   [input] t Vector of size n-1
   [output] d_ext Vector of size n+1
*/
/* SAM_LISTING_BEGIN_9*/   
VectorXd compute_d (const VectorXd &c, const VectorXd &t) {
  const unsigned int n = c.size(); //number of intervals
  VectorXd d_ext(n + 1);
  //TO DO (): compute coefficients d_j for j = 0...n 
  //Hint: periodic conditions give d_0 = d_n
  //START
  
  //END
  return d_ext;
}
/* SAM_LISTING_END_9 */
  
  
/* [input] t  Vector of size n-1
   [input] x  Vector of size N
   [input] y  Vector of size n
   [output] fval Vector of size N
*/  
/* SAM_LISTING_BEGIN_2 */
VectorXd quadspline(const VectorXd &t, const VectorXd &y, const VectorXd &x) {
  VectorXd fval(x.size());
  
  //TO DO (6-6.h): evaluate the spline at the points defined in x.
  // Assume that x is sorted.
  assert(x.minCoeff() >= 0 && x.maxCoeff() <= 1 
                && "evaluation samples out of range");
  
  //START
  
  //END
  return fval;
}
/* SAM_LISTING_END_2 */




/* [input] filename, string
   [output] plot
*/
/* SAM_LISTING_BEGIN_3 */
void plotquadspline(const std::string &filename) {
  
  VectorXd mesh = VectorXd::LinSpaced(9,.1,.9);
  auto f = [] (VectorXd t) {
    return (2 * M_PI * t).array().sin().exp();
  };
  // TO DO (6-6.i) : plot the quadratic spline for the function f based on
  // the intervals defined in t. Plot also the data (interpolation) points.
  
  //START
  
  //END
}
/* SAM_LISTING_END_3 */





/* [input] q, integer
   [output] Err vector of size n = 2^q
*/
/* SAM_LISTING_BEGIN_4 */
std::vector<double> qsp_error (unsigned int q) {
  
  std::vector<double> Err;
  unsigned int N = 10000; 
  VectorXd x = VectorXd::LinSpaced(N, 0, 1);
  auto f = [] (VectorXd t) {
    return (2 * M_PI * t).array().sin().exp();
  };
  //TO DO (6-6.j) : compute L^infty errors for all n = 2, 4, ..., 2^q
  //START 
  
  //END
  return Err;
}
/* SAM_LISTING_END_4 */



#endif
