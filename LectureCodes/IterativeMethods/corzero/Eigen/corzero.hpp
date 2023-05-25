///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>

namespace corzero {


using Eigen::VectorXd;

inline
/* SAM_LISTING_BEGIN_0 */
void fpit(double x0, VectorXd &rates, VectorXd &err)
{
  const Eigen::Index N = 15;
  double x = x0; // \com{initial guess}
  VectorXd y(N);

  for (int i=0; i<N; ++i) {
    x = x + (cos(x)+1)/sin(x);
    y(i) = x;
  }
  err.resize(N); rates.resize(N);
  err = y-VectorXd::Constant(N,x);
  rates = err.bottomRows(N-1).cwiseQuotient(err.topRows(N-1));
}
/* SAM_LISTING_END_0 */


} //namespace corzero