///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>,
///            Julien Gacon <jgacon@ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
/* SAM_LISTING_BEGIN_0 */
#include "matplotlibcpp.h" // Tools for plotting, see https://github.com/lava/matplotlib-cpp
#include <Eigen/Dense>
#include <string>
namespace plt = matplotlibcpp;
using namespace Eigen;

// Produce spy-plot of a dense Eigen matrix.
void spy(const Eigen::MatrixXd &M, const std::string &fname) {
  plt::figure();
  plt::spy(M, {{"marker", "o"}, {"markersize", "2"}, {"color", "b"}});
  plt::title("nnz = " + std::to_string(M.nonZeros()));
  plt::savefig(fname);
}

int main() {
  int n = 100;
  MatrixXd A(n, n), B(n, n);
  A.setZero();
  B.setZero();
  // Initialize matrices, see \cref{mmstruc3} 
  A.diagonal() = VectorXd::LinSpaced(n, 1, n);
  A.col(n - 1) = VectorXd::LinSpaced(n, 1, n);
  A.row(n - 1) = RowVectorXd::LinSpaced(n, 1, n);
  B = A.colwise().reverse();
  // Matrix products
  MatrixXd C = A * A, D = A * B;
  spy(A, "Aspy_cpp.eps"); // Sparse arrow matrix
  spy(B, "Bspy_cpp.eps"); // Sparse arrow matrix
  spy(C, "Cspy_cpp.eps"); // Fully populated matrix
  spy(D, "Dspy_cpp.eps"); // Sparse "framed" matrix
  return 0;
}
/* SAM_LISTING_END_0 */
