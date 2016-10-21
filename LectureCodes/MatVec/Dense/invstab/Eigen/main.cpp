///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <figure/figure.hpp>
using namespace std;
using namespace Eigen;

int main () {
/* SAM_LISTING_BEGIN_0 */
int n = 10;
VectorXd u = VectorXd::LinSpaced(n,1,n) / 3.0;
VectorXd v = u.cwiseInverse().array() * VectorXd::LinSpaced(n,1,n).unaryExpr( [](double x){return pow(-1,x);} ).array();
VectorXd x = VectorXd::Ones(n);
VectorXd expo = VectorXd::LinSpaced(19,-5,-14);
MatrixXd res(expo.size(),4);
for(int i = 0; i < expo.size(); ++i){
	double epsilon = std::pow(10, expo(i));
	MatrixXd A = u*v.transpose() + epsilon * (MatrixXd::Random(n,n) + MatrixXd::Ones(n,n))/2;
	VectorXd b = A * x;
	double nb = b.lpNorm<Infinity>();
	VectorXd xt = A.lu().solve(b);	// stable solving
	VectorXd r = b - A*xt;	// residual
	MatrixXd B = A.inverse();
	VectorXd xi = B*b;	// solving via inverse
	VectorXd ri = b - A*xi; // residual
	MatrixXd R = MatrixXd::Identity(n,n) - A*B; // residual
	res(i,0) = epsilon; res(i,1) = (r).lpNorm<Infinity>()/nb;
	res(i,2) = ri.lpNorm<Infinity>()/nb;
	// L-infinity condition number
	res(i,3) = R.lpNorm<Infinity>() / B.lpNorm<Infinity>();
}
/* SAM_LISTING_END_0 */
	// Plotting (ugly)
	mgl::Figure fig;
	fig.setlog(true, true);
	fig.plot(res.col(0),res.col(1), " +r-").label("lu() solving");
	fig.plot(res.col(0),res.col(2), " *b-").label("multiplication with inverse");
	fig.plot(res.col(0),res.col(3), " ^m").label("inverse");
	fig.xlabel("\\epsilon");// not nice
	fig.ylabel("relative residual");
	fig.grid(true);
	fig.legend();
	fig.save("invstab");
	return 0;
}
