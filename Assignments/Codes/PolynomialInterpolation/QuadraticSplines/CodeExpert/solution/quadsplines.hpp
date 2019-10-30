#ifndef QSPLINES_HPP
#define QSPLINES_HPP

//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
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
//The template will not compile until this fuction is implemented 
/* SAM_LISTING_BEGIN_0 */
std::pair<VectorXd, VectorXd> increments(const VectorXd &t) {
  // TO DO (6-6.g) : compute the increments t_j - t_{j-1} 
  // for j = 0,...,n+1 and t_{j+1} - t_{j-1} for j = 0,...,n
  // Assume that t is sorted and that t does not include the endpoints.
  // Hint: use periodic conditions to define t_{-1}, t_{n+1} for vectorization.
  
  //START
  // Create extended vector using the periodicity:
  const unsigned int n = t.size() + 1; //number of intervals
  VectorXd ext_t(n+3);
  ext_t << t(n-2) - 1, 0, t, 1, 1 + t(0) ;
  //Increments in t
  VectorXd dt = ext_t.tail(n+2).array() - ext_t.head(n+2).array();
  VectorXd ddt = ext_t.tail(n+1).array() - ext_t.head(n+1).array();
  
  return std::make_pair(dt,ddt);
  //END
}
/* SAM_LISTING_END_0 */




/* [input] t  Vector of size n-1
   [input] y  Vector of size n
   [output] c Vector of size n
*/
/* SAM_LISTING_BEGIN_1 */
VectorXd compute_c( const VectorXd &t, const VectorXd &y) {
  const unsigned int n = y.size(); //number of intervals
  VectorXd c(n);
  // TO DO (6-6.g) : Build the (sparse) matrix for spline interpolation 
  // at midpoints. Then compute the coefficients c for the data y. 
  
  //START
  assert(t.size() == n-1 && "number of intervals mismatch");
  assert(t.minCoeff() > 0 && t.maxCoeff() < 1 
                && "mesh nodes out of range");
  
  std::pair<VectorXd, VectorXd> p = increments(t);
  VectorXd dt = p.first;
  VectorXd ddt = p.second;
  
  //Note: A_0 = A_n and C_1 = C_{n+1}
  VectorXd A = dt.tail(n).cwiseQuotient(2 * ddt.tail(n));
  VectorXd C = dt.head(n).cwiseQuotient(2 * ddt.head(n));
  
  // std::vector< Triplet<double> > triplets(3*n);
  SparseMatrix<double> M(n,n);
  M.reserve( RowVectorXi::Constant(n, 3) );
  //Build matrix as triplets
  for (unsigned int j = 0; j < n - 1 ; ++j ) {
    //triplets.push_back(Triplet<double>(j, j, A(j) + C(j) + 1 ));
    //triplets.push_back(Triplet<double>(j , j + 1, C(j+1) ));
    //triplets.push_back(Triplet<double>(j + 1, j, A(j)  ));
    M.insert( j , j ) = A(j) + C(j) + 1 ;
    M.insert( j, j + 1) = C(j + 1);
    M.insert( j + 1, j ) = A(j);
  }
  // triplets.push_back(Triplet<double>(n - 1, n - 1, A(n - 1) + C(n - 1) + 1 ));
  // triplets.push_back(Triplet<double>(n - 1, 0, C(0) ));
  // triplets.push_back(Triplet<double>(0, n - 1, A(n - 1)));
  M.insert( n - 1 , n - 1 ) = A(n - 1) + C(n - 1) + 1 ;
  M.coeffRef( n - 1, 0) += C(0); // for n = 2 we edit M(1,0) and M(0,1)
  M.coeffRef( 0, n - 1 ) += A(n - 1);
  
  //M.setFromTriplets(triplets.begin(), triplets.end());
  M.makeCompressed();
  //solve linear system
  SparseLU<SparseMatrix<double>> solver;
  solver.compute(M);
  c = solver.solve(y);
  //END
  return c;
}  
/* SAM_LISTING_END_1 */
  
  
  
  
  
/* [input] c  Vector of size n
   [input] t Vector of size n-1
   [output] d_ext Vector of size n+1
*/
/* SAM_LISTING_BEGIN_9*/   
VectorXd compute_d ( const VectorXd &c, const VectorXd &t ) {
  const unsigned int n = c.size(); //number of intervals
  VectorXd d_ext(n + 1);
  //TO DO (): compute coefficients d_j for j = 0...n 
  //Hint: periodic conditions give d_0 = d_n
  //START
  std::pair<VectorXd, VectorXd> p = increments(t);
  VectorXd dt = p.first;
  VectorXd ddt = p.second;
  
  //extend c periodically for vectorization
  VectorXd c_ext(c.size() + 1);
  c_ext << c , c(0);
  
  //coefficients d_j for j = 1...n
  VectorXd d = c.cwiseProduct(dt.tail(n)) + 
                c_ext.tail(n).cwiseProduct(dt.segment(1,n));
  d = 2 * d.cwiseQuotient(ddt.tail(n));
  //extend d by periodic condition
  d_ext << d(n - 1), d; 
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
  VectorXd c = compute_c(t, y);
  const unsigned int n = y.size(); //number of intervals
  
  std::pair<VectorXd,VectorXd> p = increments(t);
  VectorXd dt = p.first;
  VectorXd ddt = p.second;
  
  VectorXd d_ext = compute_d(c, t);
  //extend t for vectorization
  VectorXd t_ext(n + 1);
  t_ext << 0, t, 1;
  
  unsigned int i = 0;
  //loop over intervals
  for (unsigned int j = 1; j <= n; ++j) {
    //loop over sample points in the interval 
    while ( i < x.size() && x(i) <= t_ext(j) ) {
      double tau = (x(i) - t_ext(j - 1)) / dt(j); 
      double uat = 1 - tau;
      fval(i) = d_ext(j) * tau * tau + 4 * c(j - 1) * tau * uat 
                      + d_ext(j - 1) * uat * uat;
      ++i;
    }
    if (i == x.size()) {break;}
  }
  //END
  return fval;
}
/* SAM_LISTING_END_2 */




/* [input] filename, string
   [output] plot
*/
/* SAM_LISTING_BEGIN_3 */
void plotquadspline( const std::string &filename ) {
  
  VectorXd mesh = VectorXd::LinSpaced(9,.1,.9);
  auto f = [] (VectorXd t) {
    return (2 * M_PI * t).array().sin().exp();
  };
  // TO DO (6-6.i) : plot the quadratic spline for the function f based on
  // the intervals defined in t. Plot also the data (interpolation) points.
  
  //START
  plt::figure();
  // The interpolation points are the midpoints of the intervals
  VectorXd t = VectorXd::LinSpaced(10, .05, .95) ;
  VectorXd y = f(t);
  // plot data points
  plt::plot(t, y, "o",{{"label","data"}} );
  
  // std::cout << y <<std::endl;
  VectorXd x = VectorXd::LinSpaced(200,0,1);
  VectorXd spline_val = quadspline(mesh, y, x);
  // std::cout << spline_val <<std::endl;
  
  plt::plot(x, spline_val, {{"label","spline"}} );
  plt::title("Quadratic spline");
  plt::xlabel("t");
  plt::ylabel("s(t)");

  plt::savefig("./cx_out/" + filename + ".png");
  //END
}
/* SAM_LISTING_END_3 */





/* [input] q, integer
   [output] Err vector of size n = 2^q
*/
/* SAM_LISTING_BEGIN_4 */
std::vector<double> qsp_error(unsigned int q) {
  
  std::vector<double> Err;
  unsigned int N = 10000; 
  VectorXd x = VectorXd::LinSpaced(N, 0, 1);
  auto f = [] (VectorXd t) {
    return (2 * M_PI * t).array().sin().exp();
  };
  //TO DO (6-6.j) : compute L^infty errors for all n = 2, 4, ..., 2^q
  //START 
  unsigned int n_max = std::pow(2,q);
  
  for(unsigned int n = 2; n <= n_max; n *= 2 ) {
    VectorXd t = VectorXd::LinSpaced(n - 1, 1.0/n, 1.0 - 1.0/n);
    // The interpolation points are the midpoints of the intervals
    VectorXd interp = VectorXd::LinSpaced(n, .5/n, 1.0 - .5/n);
    VectorXd y = f(interp);
    VectorXd f_spline = quadspline(t, y, x);
    VectorXd f_exact = f(x);
    // compute max error
    Err.push_back( (f_spline - f_exact).cwiseAbs().maxCoeff() );
  }
  //END
  return Err;
}
/* SAM_LISTING_END_4 */



#endif