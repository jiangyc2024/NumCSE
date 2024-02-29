///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

// **********************************************************************
// Eigen codes for testing sparse linear algebra routines
// **********************************************************************

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>


using std::cout;
using std::endl;
using std::cerr;
using Eigen::Triplet;
using Eigen::SparseMatrix;
using Eigen::RowMajor;
using Eigen::VectorXd;
using Eigen::SparseLU;

/* SAM_LISTING_BEGIN_0 */
template <class SpMat>
SpMat initSparseMatrix(int n) {
  using scalar_t = typename SpMat::Scalar;
  using index_t = typename SpMat::Index;
  std::vector<Triplet<scalar_t>> triplets(5 * n);

  for (index_t l = 0; l < n; ++l) {
    triplets.emplace_back(l, l, static_cast<scalar_t>(5.0));
  }
  for (index_t l = 1; l < n; ++l) {
    triplets.emplace_back(l - 1, l, static_cast<scalar_t>(1.0));
    triplets.emplace_back(l, l - 1, static_cast<scalar_t>(1.0));
  }
  const index_t m = n / 2;
  for (index_t l = 0; l < m; ++l) {
    triplets.emplace_back(l, l + m, static_cast<scalar_t>(1.0));
    triplets.emplace_back(l + m, l, static_cast<scalar_t>(1.0));
  }
  SpMat M(n, n);
  M.setFromTriplets(triplets.begin(), triplets.end());
  return M;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <class SpMat>
void printTriplets(const SpMat &M) {
  for (int k = 0; k < M.outerSize(); ++k) {
    for (typename SpMat::InnerIterator it(M, k); it; ++it){
      cout << "(" << it.row() << ',' << it.col() << ") -> " << it.value()
           << endl;
    }
  }
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
int main() {
  const int n(100);
  using SpMat = SparseMatrix<double, RowMajor>;

  const auto M = initSparseMatrix<SpMat>(n);
  cout << "M = " << M.rows() << 'x' << M.cols() << "-matrix with "
       << M.nonZeros() << " non-zeros" << endl;
  printTriplets(M);

  const VectorXd b = VectorXd::Random(n);
  VectorXd x(n);

  SparseLU<SpMat> solver;
  solver.compute(M);
  if (solver.info() != Eigen::Success) {
    cerr << "Decomposition failed!" << endl;
    return 1;
  }
  x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    cerr << "Solver failed!" << endl;
    return 1;
  }
  cout << "Residual norm = " << (M * x - b).norm() << endl;
}
/* SAM_LISTING_END_2 */
