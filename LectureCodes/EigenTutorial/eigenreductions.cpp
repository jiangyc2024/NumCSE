/* Demonstraction code for coursed Numerical Methods  for CSE, ETH Zurich
   Reductions  in Eigen, see
   https://eigen.tuxfamily.org/dox/group__TutorialReductionsVisitorsBroadcasting.html
   @author Ralf Hiptmair
   @date August 2020
*/

#include <Eigen/Dense>
#include <iostream>

/* SAM_LISTING_BEGIN_R */
template <class Matrix> void sumEntries(Eigen::MatrixBase<Matrix> &M) {
  using Scalar = typename Eigen::MatrixBase<Matrix>::Scalar;
  // Compute sum  of all entries
  const Scalar s = M.sum();
  // Row-wise and column-wise sum of entries: results are vectors
  Eigen::Matrix<Scalar, 1, Eigen::Dynamic> colsums{M.colwise().sum()};
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> rowsums{M.colwise().sum()};
  std::cout << M.rows() << 'x' << M.cols() << "-matrix: " << colsums.sum()
            << " = " << rowsums.sum() << " = " << s << std::endl;
}
/* SAM_LISTING_END_R */

int main(int /*argc*/, char ** /*argv*/) {
  sumEntries((Eigen::Matrix3d() << 1, 2, 3, 4, 5, 6, 7, 8, 9).finished());
  sumEntries(Eigen::MatrixXd(3, 4).setOnes());
  return 0;
}

