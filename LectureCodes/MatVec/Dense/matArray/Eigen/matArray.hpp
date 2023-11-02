///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>
#include <iostream>

namespace matarray {


using std::cout;
using std::endl;

/* SAM_LISTING_BEGIN_0 */
void inline matArray(int nrows,int ncols){
  Eigen::MatrixXd m1(nrows,ncols);
  Eigen::MatrixXd m2(nrows,ncols);
  for(int i = 0; i < m1.rows(); i++) {
    for(int j = 0; j < m1.cols(); j++) {  
      m1(i,j) = static_cast<double>((i+1))/(j+1); 
      m2(i,j) = static_cast<double>((j+1))/(i+1); 
    }
  }
  // \com{Entry-wise} product, not a matrix product
  const Eigen::MatrixXd m3 = (m1.array() * m2.array()).matrix();
  // Explicit entry-wise operations on matrices are possible
  const Eigen::MatrixXd m4(m1.cwiseProduct(m2)); 
  // \com{Entry-wise} logarithm 
  cout << "Log(m1) = " << endl << log(m1.array()) << endl;
  // Entry-wise boolean expression, \texttt{true} cases counted
  cout << (m1.array() > 3).count() << " entries of m1 > 3" << endl;
}
/* SAM_LISTING_END_0 */


} //namespace matarray