///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
template<typename MatType>
void reshapetest(MatType &M)
{
  using index_t = typename MatType::Index;
  using entry_t = typename MatType::Scalar;
  const index_t nsize(M.size());

 // reshaping possible only for matrices with non-prime dimensions 
  if ((nsize %2) == 0) {
    entry_t *Mdat = M.data(); // raw data array for M
    // Reinterpretation of data of M
    Map<Eigen::Matrix<entry_t,Dynamic,Dynamic>> R(Mdat,2,nsize/2);
    // (Deep) copy data of M into matrix of different size
    Eigen::Matrix<entry_t,Dynamic,Dynamic> S =
      Map<Eigen::Matrix<entry_t,Dynamic,Dynamic>>(Mdat,2,nsize/2);

    cout << "Matrix M = " << endl << M << endl;
    cout << "reshaped to " << R.rows() << 'x' << R.cols()
	 <<  " = " << endl << R << endl;
    // Modifying R affects M, because they share the data space !
    R *= -1.5;
    cout << "Scaled (!) matrix M = " << endl << M << endl;
    // Matrix S is not affected, because of deep copy
    cout << "Matrix S = " << endl << S << endl;
   }
}
/* SAM_LISTING_END_0 */
