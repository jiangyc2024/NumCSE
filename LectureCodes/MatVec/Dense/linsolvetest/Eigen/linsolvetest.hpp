#pragma once

#include <vector>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! Eigen code: assessing the gain from using special properties
//! of system matrices in Eigen
MatrixXd timing(){
	std::vector<int> n = {16,32,64,128,256,512,1024,2048,4096,8192};
	int nruns = 3;
	MatrixXd times(n.size(),3);
	for(int i = 0; i < n.size(); ++i){
		Timer t1, t2;	// timer class
		MatrixXd A = VectorXd::LinSpaced(n[i],1,n[i]).asDiagonal() ;
		A += MatrixXd::Ones(n[i],n[i]).triangularView<Upper>();
		VectorXd b = VectorXd::Random(n[i]);
		VectorXd x1(n[i]), x2(n[i]);
		for(int j = 0; j < nruns; ++j){
			t1.start();	x1 = A.lu().solve(b);	t1.stop();
			t2.start();	x2 = A.triangularView<Upper>().solve(b); t2.stop();
		}
		times(i,0) = n[i]; times(i,1) = t1.min(); times(i,2) = t2.min();
	}
	return times;
}
/* SAM_LISTING_END_0 */
