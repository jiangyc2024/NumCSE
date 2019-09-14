///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;

int main() {
  MatrixXd A(3, 3);
  A.setRandom();
  /* SAM_LISTING_BEGIN_0 */
  const Eigen::MatrixXd::Index n = A.cols();
  assert(n == A.rows()); // ensure square matrix
  Eigen::PartialPivLU<MatrixXd> lu(A);
  // Normalized lower-triangule factor
  MatrixXd L = MatrixXd::Identity(n, n);
  L.triangularView<StrictlyLower>() += lu.matrixLU();
  // Upper triangular factor
  MatrixXd U = lu.matrixLU().triangularView<Upper>();
  // Permutation matrix, see \cref{def:permmat}
  MatrixXd P = lu.permutationP();
  /* SAM_LISTING_END_0 */
  std::cout << "A = " << A << std::endl;
  std::cout << "L = " << L << std::endl;
  std::cout << "U = " << U << std::endl;
  // Verify correctness of decomposition
  std::cout << "residual matrix = " << P * L * U - A << std::endl;
  return 0;
}
