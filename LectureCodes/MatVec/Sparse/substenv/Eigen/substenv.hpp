#pragma once
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#pragma begin<0>
//! evelope aware forward substitution for \Blue{$\VL\Vx=\Vy$}
//! (\Blue{$\VL$} = lower triangular matrix)
//! argument \texttt{mc}: column bandwidth vector------------------------------------
VectorXd substenv(const MatrixXd &L, const VectorXd &y, const VectorXi &mr){
	int n = L.cols(); VectorXd x(n);
	x(0) = y(0)/L(0,0);
	for(int i = 1; i < n; ++i){
		if( mr(i) > 0){
			double zeta = L.row(i).segment(i-mr(i),mr(i))
					* x.segment(i-mr(i), mr(i));
			x(i) = (y(i) - zeta)/L(i,i);
		}
		else x(i) = y(i) / L(i,i);
	}
	return x;
}
#pragma end<0>
