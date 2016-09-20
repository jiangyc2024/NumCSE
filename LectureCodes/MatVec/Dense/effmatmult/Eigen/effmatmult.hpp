///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <iomanip>
#include <cmath>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! This function compares the runtimes for the multiplication
//! of a vector with a rank-1 matrix $\Va\Vb^{\top}$, $\Va,\Vb\in\bbR^n$
//! using different associative evaluations.
//! Runtime measurements consider minimal time for
//! several (\texttt{nruns}) runs
MatrixXd dottenstiming() {
  const int nruns = 3, minExp = 2, maxExp = 13;
  // Matrix for storing recorded runtimes
  MatrixXd timings(maxExp-minExp+1,3);	
  for(int i = 0; i <= maxExp-minExp; ++i){
    Timer tfool, tsmart; // Timer objects
    const int n = std::pow(2, minExp + i);
    VectorXd a = VectorXd::LinSpaced(n,1,n);
    VectorXd b = VectorXd::LinSpaced(n,1,n).reverse();
    VectorXd x = VectorXd::Random(n,1), y(n);
    for(int j = 0; j < nruns; ++j){
      // Grossly wasteful evaluation
      tfool.start(); y = (a*b.transpose())*x;	tfool.stop();
      // Efficient implementation
      tsmart.start(); y = a * b.dot(x);	tsmart.stop();	      
    }
    timings(i,0)=n;
    timings(i,1)=tsmart.min(); timings(i,2)=tfool.min();
  }
  return timings;
}
/* SAM_LISTING_END_0 */
