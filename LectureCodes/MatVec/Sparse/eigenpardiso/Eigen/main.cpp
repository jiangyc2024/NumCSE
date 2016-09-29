///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

using namespace std;

/* SAM_LISTING_BEGIN_0 */
template <class SpMat>
SpMat initSparseMatrix(size_t n)
{
  using index_t = typename SpMat::Index;
  using scalar_t = typename SpMat::Scalar;
  vector<Eigen::Triplet<scalar_t>> triplets(5*n);

  for (size_t l=0; l<n; ++l)
    triplets.push_back(Eigen::Triplet<scalar_t>(l,l,5.0));
  for (size_t l=1; l<n; ++l) {
    triplets.push_back(Eigen::Triplet<scalar_t>(l-1,l,1.0));
    triplets.push_back(Eigen::Triplet<scalar_t>(l,l-1,1.0));
  }
  const size_t m = n/2;
  for (size_t l=0; l<m; ++l) {
    triplets.push_back(Eigen::Triplet<scalar_t>(l,l+m,1.0));
    triplets.push_back(Eigen::Triplet<scalar_t>(l+m,l,1.0));
  }
  SpMat M(n,n);
  M.setFromTriplets(triplets.begin(), triplets.end());
  M.makeCompressed();
  return M;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_2 */
void solveSparsePardiso(size_t n){
  using SpMat = Eigen::SparseMatrix<double>;
  // Initialize a sparse matrix
  const SpMat M = initSparseMatrix<SpMat>(n);
  const Eigen::VectorXd b = Eigen::VectorXd::Random(n);
  Eigen::VectorXd x(n);
  // Initalization of the sparse direct solver based on the Pardiso library with directly passing the matrix M to the solver
  // Pardiso is part of the Intel MKL library, see also \cref{ex:mmeigenmkl}
  Eigen::PardisoLU<SpMat> solver(M); 
  // The checks of \cref{cpp:splsetest} are omitted
  // solve the LSE
  x = solver.solve(b);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
int main(){
  solveSparsePardiso(100);
  return 0;
}
/* SAM_LISTING_END_3 */
