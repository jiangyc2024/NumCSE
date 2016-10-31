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
// minimizing the mean square error of Euclidean distances. The columns
// of the matrix X contain the point coordinates, n specifies the 
// desired number of clusters.
std::pair<MatrixXd, VectorXi> pointcluster(const MatrixXd & X, const int n){
  int N = X.cols(); // no. of points
  int k = X.rows(); // dimension of space
  
  // Start with two clusters obtained by principal axis separation
  int nc = 1; // Current number of clusters
  // Initial single cluser encompassing all points
  VectorXi Ibig = VectorXi::LinSpaced(N,1,N);
  int nbig = 0; // Index of largest cluster
  MatrixXd C = X.rowwise().sum()/N; // center of gravity
  VectorXi idx(N); idx.setOnes();
  // Split largest cluster into two using the principal axis separation
  // algorithm
  while(nc < n){
	  VectorXi i1, i2;
	  MatrixXd Xbig(k,Ibig.size());
	  for(int i = 0; i < Ibig.size(); ++i)	// slicing
		Xbig.col(i) = X.col(Ibig(i));
	  std::tie(i1, i2) = princaxissep(Xbig);
	  
	  VectorXi i1_tmp(i1.size()), i2_tmp(i2.size());
	  VectorXd c1(k), c2(k);
	  c1.setZero(); c2.setZero();
	  for(int i = 0; i < i1.size(); ++i){
		  i1_tmp(i) = Ibig(i1(i));/// CAN BE SIMPLIFIED
		  c1 += X.col(i1_tmp(i));
	  }
	  i1 = i1_tmp;
	  for(int i = 0; i < i2.size(); ++i){
		  i2_tmp(i) = Ibig(i2(i));
		  c2 += X.col(i2_tmp(i));
	  }
	  c1 /= i1.size();
	  c2 /= i2.size();
	  assert(i1.size() >= i2.size()); /// does this condition always hold?
	  //std::cout << "Test1" << std::endl;
	  C.col(nbig) = c1;
	  C.conservativeResize(k,nc+1); /// might not be necessary if passing only blocks of array
	  //std::cout << "Test2" << std::endl;
	  C.rightCols<1>() = c2;
	  ++nc;
	  // Improve clusters by Lloyd-Max iteration
	  //VectorXi idx;
	  VectorXd cds;
	  //std::cout << X << std::endl << C << std::endl;
	  
	  //std::cout << "before:\n" << idx << std::endl;
	  lloydmax(X,C,idx,cds);
	  std::cout << "after:\n" << idx << std::endl;
	  // Identify cluster with biggest contribution to mean square error
	  double cdm = cds.maxCoeff(&nbig);
	  std::cout << "cds\n" << cds << std::endl;
	  int counter = 0;
	  for(int i = 0; i < idx.size(); ++i){
		  if(idx(i) == nbig){
			  Ibig(counter) = i;
			  ++counter;
		  }
	  }
	  Ibig.conservativeResize(counter);
	  std::cout << "Ibig\n" << Ibig << " end"<< std::endl;
	  
  }
  
  
  return std::make_pair(C, idx);
}
/* SAM_LISTING_END_0 */



/// algorithms are writen for matlab slicing style, might be better to rewrite everything into c++ style?
