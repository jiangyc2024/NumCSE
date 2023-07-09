///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>
#include <vector>

namespace princaxissep {

using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;

inline
/* SAM_LISTING_BEGIN_0 */
// Separation of a set of points whose coordinates are stored in the 
// columns of X according to their location w.r.t. the principal axis
std::pair<VectorXi, VectorXi> princaxissep(const MatrixXd & X){
  const Eigen::Index N = X.cols();	// no. of points
  const VectorXd g = X.rowwise().sum() / N; // Center of gravity, \emph{cf.} \eqref{eq:cgrav}
  MatrixXd Y = X - g.replicate(1,N); // Normalize point coordinates.
  // Compute \Red{principal axes}, \emph{cf.} \eqref{eq:pax} and \eqref{lsq:maxconst}. Note that the
  // SVD of a symmetric matrix is available through an orthonormal 
  // basis of eigenvectors.
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(Y*Y.transpose());
  // Major principal axis
  Eigen::VectorXd a = es.eigenvectors().rightCols<1>();
  // Coordinates of points w.r.t. to major principal axis
  Eigen::VectorXd c = a.transpose()*Y;
  // Split point set according to locations of projections on principal axis
  // std::vector with indices to prevent resizing of matrices
  std::vector<int> i1;
  std::vector<int> i2;
  for(int i = 0; i < c.size(); ++i){
    if(c(i) >= 0) {
      i1.push_back(i);
    }
    else {
      i2.push_back(i);
    }
  }
  // return the mapped std::vector as Eigen::VectorXd
  return {
    VectorXi::Map(i1.data(), static_cast<Eigen::Index>(i1.size())), 
    VectorXi::Map(i2.data(), static_cast<Eigen::Index>(i2.size()))
  };
}
/* SAM_LISTING_END_0 */


} //namespace princaxissep