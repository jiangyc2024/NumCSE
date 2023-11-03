///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "spy.hpp"
#include <Eigen/Sparse>
#include <iostream>
#include <span>
#include <string>
#include <unsupported/Eigen/SparseExtra>  // for import

using Eigen::SparseMatrix;
using Eigen::SimplicialLDLT;
using Eigen::Lower;
using Eigen::AMDOrdering;
using Eigen::NaturalOrdering;
using Eigen::MatrixXd;

//NOLINTBEGIN(bugprone-exception-escape)
int main(int argc, char *argv[]) {
  const std::span args(argv, argc);
  typedef SparseMatrix<double> SpMat_t;
  SpMat_t M;
  // load file
  std::string filename = "poisson2D.mtx";
  if (argc >= 2) {
    filename = args[1];
  }
  if (! loadMarket(M, "./" + filename)) {
    std::cout << "failed import of " << filename << std::endl;
    return 0;
  }
  std::cout << "import successful\n";
  /* SAM_LISTING_BEGIN_0 */
  // L and U cannot be extracted from SparseLU --> LDLT
  const SimplicialLDLT<SpMat_t, Lower, AMDOrdering<int> > solver1(M);
  const SimplicialLDLT<SpMat_t, Lower, NaturalOrdering<int> > solver2(M);
  const MatrixXd U1 =
      solver1.matrixU() *
      MatrixXd::Identity(
          M.rows(), M.cols());  // explicit conversion fixes occasional segfault
  const MatrixXd U2 = MatrixXd(solver2.matrixU());
  // Plotting
  spy(M, "Sparse matrix M", "MSpy.eps");
  spy(U1, "U factor (approximate minimum degree)", "AMDUSpy.eps");
  spy(U2, "U factor (no reordering)", "NaturalUSpy.eps");
  /* SAM_LISTING_END_0 */
  return 0;
}
//NOLINTEND(bugprone-exception-escape)
