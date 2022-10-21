///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>
#include <cmath>

namespace ApproxPIinstable {


using std::sqrt;
using Eigen::MatrixXd;

inline
/* SAM_LISTING_BEGIN_0 */
//! Approximation of Pi by approximating the circumference of a
//! regular polygon
MatrixXd ApproxPIinstable(double tol = 1e-8, unsigned int maxIt = 50){
  double s=sqrt(3)/2.; 
  double An=3.*s;// initialization (hexagon case)
  unsigned int n = 6;
  unsigned int it = 0;				
  MatrixXd res(maxIt,4);	// matrix for storing results
  res(it,0) = n; res(it,1) = An;
  res(it,2) = An - M_PI; res(it,3)=s;
  while( it < maxIt && s > tol ){// terminate when s is ’small enough’
  	s = sqrt((1.- sqrt(1.-s*s))/2.);// recursion for area \Label{pis:1}
  	n *= 2; An = n/2.*s;	// new estimate for circumference
  	++it;
  	res(it,0) =n; res(it,1) =An;// store results and (absolute) error
  	res(it,2) = An - M_PI; res(it,3)=s; 	
  }
  return res.topRows(it);
}
/* SAM_LISTING_END_0 */


} //namespace ApproxPIinstable
