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
std::tuple<double, VectorXi, VectorXd> distcomp(const MatrixXd & X, const MatrixXd & C){
  // Compute squared distances
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
	  mx(j) = d.col(j).minCoeff(&idx(j));///
  }
  double sumd = mx.sum();
  // Computer sum of squared distances within each cluster
  VectorXd cds(C.cols()); cds.setZero();
  for(int j = 0; j < idx.size(); ++j)
	  cds(idx(j)) += mx(j);
  std::cout << "distcomp \n" << idx << "  end" << std::endl;
  return std::make_tuple(sumd, idx, cds);
}

// Lloyd-Max iterative vector quantization algorithm for discrete point 
// sets; the columns of \texttt{X} contain the points \Blue{$\Vx_{i}$}, the columns of 
// \texttt{C} initial approximations for the centers of the clusters. The final 
// centers are returned in \texttt{C}, the index vector \texttt{idx} specifies 
// the association of points with centers.
void lloydmax(const MatrixXd & X, MatrixXd & C, VectorXi & idx, VectorXd & cds, const double tol = 0.0001){
  int k = X.rows(); // dimension of space
  int N = X.cols();	// no. of points
  int n = C.cols(); // no. of clusters
  
  if(k != C.rows())
	throw std::logic_error("dimension mismatch");
  
  double sd_old = std::numeric_limits<double>::max();
  double sd = 1;
  std::tie(sd, idx, cds) = distcomp(X,C);
  // Terminate, if sum of squared minimal distances has not changed much
  while( (sd_old-sd)/sd > tol ){
	  // Compute new centers of gravity according to \eqref{eq:cgrav}
	  /*
	  C.setZero();
	  VectorXi nj(n); // number of points in cluster for normalization
	  for(int j = 0; j < N; ++j){	// loop over all points
		  C.col(idx(j)) += X.col(j);
		  ++nj(idx(j));
	  }
	  // Normalization
	  C.applyOnTheRight(nj.asDiagonal().inverse());	// \eigen{} optimized
	  * could also use C.rowwise().array() /= nj.array();
	  * 
	  /// does not work i one cluster gets empty!!!!!, since C gets reset and nj can't be zero
	  /// --->
	  * 
	  */
	  for(int j = 0; j < n; ++j){ // loop over all clusters
		  int nj = 0;
		  VectorXd gj(k); gj.setZero(); // temporary for CoG
		  for(int i = 0; i < N; ++i){ // loop over all points
			  if(idx(i) == j){
				  ++nj;
				  gj += X.col(i);
			  }
		  }
		  if(nj > 0)
			C.col(j) = gj/nj;
	  }
	  sd_old = sd;
	  std::tie(sd, idx, cds) = distcomp(X,C);
  }
}
/* SAM_LISTING_END_0 */



///------> generalize principal axis separation to N-Dimensions
