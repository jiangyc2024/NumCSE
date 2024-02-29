///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

using Eigen::RowVectorXd;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::StrictlyLower;
using Eigen::Upper;
using Eigen::VectorXi;

#include "bandwidth.hpp"
#include "spdiags.hpp"
#include "substenv.hpp"

int main() {
  // Build matrix
  const Eigen::Index n = 5;
  RowVectorXd diag_el(5);
  diag_el << -1, -1, 3, -1, -1;
  VectorXi diag_no(5);
  diag_no << -n, -1, 0, 1, n;
  MatrixXd B = diag_el.replicate(2 * n, 1);
  B(n - 1, 1) = 0;
  B(n, 3) = 0;  // delete elements
  const SparseMatrix<double> A = spdiags(B, diag_no, 2 * n, 2 * n);
  // It's not possible to get L, U with the sparseLU --> dense
  auto solver = MatrixXd(A).lu();
  MatrixXd L = MatrixXd::Identity(2 * n, 2 * n);
  L += solver.matrixLU().triangularView<StrictlyLower>();
  const MatrixXd U = solver.matrixLU().triangularView<Upper>();

  std::cout << "L\n" << L << std::endl;
  const VectorXd y = VectorXd::LinSpaced(2 * n, 1, 2 * n);
  std::cout << "y\n" << y << std::endl;
  const VectorXd x_ex = L.lu().solve(y);
  std::cout << "x from Eigen\n" << x_ex << std::endl;
  const VectorXi mr = bandwidth::rowbandwidth(L);
  const VectorXd x_own = envelope::substenv(L, y, mr);
  std::cout << "x from own impl.\n" << x_own << std::endl;
  return 0;
}
