#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <figure/figure.hpp>
using namespace std;
using namespace Eigen;

int main () {
	/* SAM_LISTING_BEGIN_0 */
	MatrixXd res(100,2);
	for(int n = 10; n <= 100*10; n += 10){
		MatrixXd A(n,n); A.setIdentity();
		A.triangularView<StrictlyLower>().setConstant(-1);
		A.rightCols<1>().setOnes();
		VectorXd x = VectorXd::Constant(n,-1).binaryExpr( VectorXd::LinSpaced(n,1,n), [](double x, double y){return pow(x,y);});
		double relerr = (A.lu().solve(A*x)-x).norm()/x.norm();
		res(n/10-1,0) = n; res(n/10-1,1) = relerr;
	}
	// ... different solver(e.g. colPivHouseholderQr()), plotting
	/* SAM_LISTING_END_0 */
	return 0;
}
