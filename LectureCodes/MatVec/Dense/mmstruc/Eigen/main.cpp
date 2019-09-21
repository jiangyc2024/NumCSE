///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>,
///            Julien Gacon <jgacon@ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
/* SAM_LISTING_BEGIN_0 */
#include <string>
#include <Eigen/Dense>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace Eigen;

void spy(const Eigen::MatrixXd& M, const std::string& fname) {
  plt::figure();
  plt::spy(M, {{"marker", "o"}, {"markersize", "2"}, {"color", "b"}});
  plt::title("nnz = " + std::to_string(M.nonZeros()));
  plt::savefig(fname);
}

int main () {
	int n = 100;
	MatrixXd A(n,n), B(n,n); A.setZero(); B.setZero();
	A.diagonal() = VectorXd::LinSpaced(n,1,n);
	A.col(n-1) = VectorXd::LinSpaced(n,1,n);
	A.row(n-1) = RowVectorXd::LinSpaced(n,1,n);
	B = A.colwise().reverse();
	MatrixXd C = A*A, D = A*B;
  spy(A, "Aspy_cpp.eps");
  spy(B, "Bspy_cpp.eps");
  spy(C, "Cspy_cpp.eps");
  spy(D, "Dspy_cpp.eps");
	return 0;
}
/* SAM_LISTING_END_0 */
