///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): R. Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Sparse>
#include <iostream>
#include <vector>

// Conversion of a sparse matrix (in CRS/CCS format) to COO/triplet format
// see
// https://stackoverflow.com/questions/28685877/convert-an-eigen-matrix-to-triplet-form-c
// http://eigen.tuxfamily.org/dox/group__TutorialSparse.html
/* SAM_LISTING_BEGIN_0 */
template <typename Scalar>
std::vector<Eigen::Triplet<Scalar>> convertToTriplets(
    Eigen::SparseMatrix<Scalar> &A) {
  // vector of triplets to be grown in the following loops
  std::vector<Eigen::Triplet<Scalar>> triplets{};
  // Loop over row/columns (depending on column/row major format
  for (int k = 0; k < A.outerSize(); ++k) {
    // Loop over inner dimension and obtain triplets corresponding
    // to non-zero entries.
    for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(A, k); it;
         ++it) {
      // Retrieve triplet data from iterator
      triplets.emplace_back(it.row(), it.col(), it.value());
    }
  }
  return triplets;
}
/* SAM_LISTING_END_0 */
