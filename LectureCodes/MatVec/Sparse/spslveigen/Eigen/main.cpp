///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

// **********************************************************************
// Eigen codes for testing dense linear algebra routines
// **********************************************************************

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;

/* SAM_LISTING_BEGIN_0 */
template <class SpMat>
SpMat initSparseMatrix(size_t n)
{
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

/* SAM_LISTING_BEGIN_1 */
template <class SpMat>
void printTriplets(const SpMat &M)
{
  for (int k=0; k<M.outerSize(); ++k)
    for (typename SpMat::InnerIterator it(M,k); it; ++it)
      cout << "(" << it.row() << ',' << it.col() << ") -> " << it.value() << endl;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void solveSparseTest(size_t n)
{
  using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  
  const SpMat M = initSparseMatrix<SpMat>(n);
  cout << "M = " << M.rows() << 'x' << M.cols() << "-matrix with " << M.nonZeros()
       << "non-zeros" << endl;
  printTriplets(M); 

  const Eigen::VectorXd b = Eigen::VectorXd::Random(n);
  Eigen::VectorXd x(n);

  Eigen::SparseLU<SpMat> solver; solver.compute(M);
  if(solver.info() != Eigen::Success) {
    cerr << "Decomposition failed!" << endl;
    return;
  }
  x = solver.solve(b);
  if(solver.info() != Eigen::Success) {
    cerr << "Solver failed!" << endl;
    return;
  }
  cout << "Residual norm = " << (M*x-b).norm() << endl;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
int main(int argc,char **argv)
{
  cout << "EIGEN DENSE LINEAR ALGEBRA CODES" << endl;
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <selection>" << endl;
    return(-1L);
  }
  else {
    const int sel = atoi(argv[1]);
    switch (sel) {
    case 1: { solveSparseTest(100); break; }
    default: { cerr << "Invalid selection" << endl; exit(-1L); }
    }
  }
  return 0;
}
/* SAM_LISTING_END_3 */
