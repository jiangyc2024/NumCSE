///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <limits>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
template <class Derived>
std::tuple<double, VectorXi, VectorXd> distcomp(const MatrixXd & X, const MatrixBase<Derived> & C){
  // Compute squared distances
  // d.row(j) = squared distances from all points in X to cluster j
  MatrixXd d(C.cols(), X.cols());
  for(int j = 0; j < C.cols(); ++j){
	  MatrixXd Dv = X - C.col(j).replicate(1, X.cols());
	  d.row(j) = Dv.array().square().colwise().sum();
  }
  // Compute minimum distance point association and sum of minimal squared distances
  VectorXi idx(d.cols());  VectorXd mx(d.cols());
  for(int j = 0; j < d.cols(); ++j){
	  // mx(j) tells the minimal squared distance of point j to the nearest cluster
	  // idx(j) tells to which cluster point j belongs
	  mx(j) = d.col(j).minCoeff(&idx(j));
  }
  double sumd = mx.sum();	// sum of all squared distances
  // Computer sum of squared distances within each cluster
  VectorXd cds(C.cols()); cds.setZero();
  for(int j = 0; j < idx.size(); ++j)	// loop over all points
	  cds(idx(j)) += mx(j);
  return std::make_tuple(sumd, idx, cds);
}

// Lloyd-Max iterative vector quantization algorithm for discrete point 
// sets; the columns of \texttt{X} contain the points \Blue{$\Vx_{i}$}, the columns of 
// \texttt{C} initial approximations for the centers of the clusters. The final 
// centers are returned in \texttt{C}, the index vector \texttt{idx} specifies 
// the association of points with centers.
template <class Derived>
void lloydmax(const MatrixXd & X, MatrixBase<Derived> & C, VectorXi & idx, VectorXd & cds, const double tol = 0.0001){
  int k = X.rows(); // dimension of space
  int N = X.cols();	// no. of points
  int n = C.cols(); // no. of clusters
  if(k != C.rows())
	throw std::logic_error("dimension mismatch");
  double sd_old = std::numeric_limits<double>::max();
  double sd;
  std::tie(sd, idx, cds) = distcomp(X,C);
  // Terminate, if sum of squared minimal distances has not changed much
  while( (sd_old-sd)/sd > tol ){
	  // Compute new centers of gravity according to \eqref{eq:cgrav}
	  MatrixXd Ctmp(C.rows(),C.cols());  Ctmp.setZero();
	  // number of points in cluster for normalization
	  VectorXi nj(n); nj.setZero();
	  for(int j = 0; j < N; ++j){	// loop over all points
		  Ctmp.col(idx(j)) += X.col(j);
		  ++nj(idx(j));	// count associated points for normalization
	  }
	  for(int i = 0; i < Ctmp.cols(); ++i){
		  if(nj(i) > 0)
			C.col(i) = Ctmp.col(i)/nj(i);	// normalization
	  }
	  sd_old = sd;
	  // Get new minimum association of the points to cluster points
	  // for next iteration
	  std::tie(sd, idx, cds) = distcomp(X,C);
  }
}
// Note: this function is needed to allow a call with an rvalue
// \&\& stands for an rvalue reference and allows rvalue arguments
// such as C.leftCols(nc) to be passed by reference (\cpp 11 feature)
template <class Derived>
void lloydmax(const MatrixXd & X, MatrixBase<Derived> && C, VectorXi & idx, VectorXd & cds, const double tol = 0.0001){
	lloydmax(X, C, idx, cds);
}
/* SAM_LISTING_END_0 */
