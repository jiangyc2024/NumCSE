///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>

#include "pas.hpp"
#include "lloydmax.hpp"

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
// n-quantization of point set in k-dimensional space based on 
// minimizing the mean square error of Euclidean distances. The 
// columns of the matrix X contain the point coordinates, n specifies 
// the desired number of clusters.
std::pair<MatrixXd, VectorXi> pointcluster(const MatrixXd & X, const int n){
  int N = X.cols(); // no. of points
  int k = X.rows(); // dimension of space
  // Start with two clusters obtained by principal axis separation
  int nc = 1; // Current number of clusters
  // Initial single cluster encompassing all points
  VectorXi Ibig = VectorXi::LinSpaced(N,0,N-1);
  int nbig = 0; // Index of largest cluster
  MatrixXd C(X.rows(), n);	// matrix for cluster midpoints
  C.col(0) = X.rowwise().sum()/N; // center of gravity
  VectorXi idx(N); idx.setOnes();
  // Split largest cluster into two using the principal axis separation
  // algorithm
  while(nc < n){
	  VectorXi i1, i2;
	  MatrixXd Xbig(k,Ibig.size());
	  for(int i = 0; i < Ibig.size(); ++i)	// slicing
		Xbig.col(i) = X.col(Ibig(i));
	  // separete Xbig into two clusters, i1 and i2 are index vectors
	  std::tie(i1, i2) = princaxissep(Xbig);
	  // new cluster centers of gravity
	  VectorXd c1(k), c2(k); c1.setZero(); c2.setZero();
	  for(int i = 0; i < i1.size(); ++i)
		  c1 += X.col(Ibig(i1(i)));
	  for(int i = 0; i < i2.size(); ++i)
		  c2 += X.col(Ibig(i2(i)));
	  c1 /= i1.size();	  c2 /= i2.size(); // normalization
	  C.col(nbig) = c1;
	  C.col(nbig+1) = c2;
	  ++nc; // Increase number of clusters
	  // Improve clusters by Lloyd-Max iteration
	  VectorXd cds;	// saves mean square error of clusters
	  // Note C.leftCols(nc) is passed as rvalue reference (\cpp 11)
	  lloydmax(X,C.leftCols(nc),idx,cds); 
	  // Identify cluster with biggest contribution to mean square error
	  cds.maxCoeff(&nbig);
	  int counter = 0;
	  // update Ibig with indices of points in cluster with biggest contribution
	  for(int i = 0; i < idx.size(); ++i){
		  if(idx(i) == nbig){
			  Ibig(counter) = i;
			  ++counter;
		  }
	  }
	  Ibig.conservativeResize(counter);
  }
  return std::make_pair(C, idx);
}
/* SAM_LISTING_END_0 */



/// algorithms are writen for matlab slicing style, might be better to rewrite everything into c++ style?
