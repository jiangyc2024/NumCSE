///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"
#include "gausselimsolve.hpp"


using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! Eigen code for timing numerical solution of linear systems
MatrixXd gausstiming(){
  std::vector<int> n = {8,16,32,64,128,256,512,1024,2048,4096,8192};
  int nruns = 3;
  MatrixXd times(n.size(),3);
  for(int i = 0; i < n.size(); ++i){
  	Timer t1, t2;	// timer class
  	MatrixXd A = MatrixXd::Random(n[i],n[i]) + n[i]*MatrixXd::Identity(n[i],n[i]);
  	VectorXd b = VectorXd::Random(n[i]);
  	VectorXd x(n[i]);
  	for(int j = 0; j < nruns; ++j){
  	  t1.start();  x = A.lu().solve(b);  t1.stop();// Eigen implementation
  	  #ifndef EIGEN_USE_MKL_ALL	// only test own algorithm without MKL
  	  if(n[i] <= 4096)		// Prevent long runs
  	  	t2.start();	gausselimsolve(A,b,x);	t2.stop();	// own gauss elimination
  	  #endif
  	}
  	times(i,0) = n[i]; times(i,1) = t1.min(); times(i,2) = t2.min();
  }
  return times;
}
/* SAM_LISTING_END_0 */
