// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::Upper;
using Eigen::Lower;

using std::cout;
using std::endl;
using std::min;

/** @brief Access to a matrix block
    @tparam Matrix any Eigen matrix compatible type
    @param M reference to a matrix-like object.

    Since Eigen uses expression templates, the expressions often evaluate to
    temporary objects, instead of Eigen::Matrix. This function demonstrates a
   way to accept such objects as arguments of functions. They are still
   compatible with Eigen::Matrix since Eigen::Matrix inherits from
   Eigen::MatrixBase.
*/
/* SAM_LISTING_BEGIN_4 */
template <typename MatType>
void blockAccess(Eigen::MatrixBase<MatType> &M) {
  using index_t = typename Eigen::MatrixBase<MatType>::Index;
  using entry_t = typename Eigen::MatrixBase<MatType>::Scalar;
  const index_t nrows(M.rows());  // No. of rows
  const index_t ncols(M.cols());  // No. of columns

  cout << "Matrix M = " << endl << M << endl << endl;  // Print matrix
  // Block size half the size of the matrix
  int p = nrows / 2;
  int q = ncols / 2;
  // Output submatrix with left upper entry at position \texttt{(i,i)}
  for (index_t i = 0; i < min(p, q); i++) {
    cout << "Block (" << i << ',' << i << ',' << p << ',' << q << ") = " << endl
         << M.block(i, i, p, q) << endl
         << endl;
  }
  // l-value access: Modify sub-matrix by adding a constant
  M.block(1, 1, p, q) += MatrixXd::Constant(p, q, static_cast<entry_t>(1));
  cout << "With modified block (1,1," << p << ',' << q << "): M = " << endl
       << M << endl
       << endl;
  // r-value access: Extract sub-matrix
  MatrixXd B(M.block(1, 1, p, q));
  cout << "Isolated modified block = " << endl << B << endl << endl;
  // Special sub-matrices
  cout << p << " top rows of M = " << endl << M.topRows(p) << endl << endl;
  cout << p << " bottom rows of M = " << endl
       << M.bottomRows(p) << endl
       << endl;
  cout << q << " left cols of M = " << endl << M.leftCols(q) << endl << endl;
  cout << q << " right cols of M = " << endl << M.rightCols(p) << endl << endl;

  // Note: for vectors only, you can also use the functions head and tail

  // r-value access to upper triangular part
  // See https://stackoverflow.com/a/24101297 as to why the ´template´ keyword
  // is required here
  const MatrixXd T = M.template triangularView<Upper>();
  cout << "Upper triangular part = " << endl << T << endl << endl;

  // l-value access to upper triangular part
  M.template triangularView<Lower>() *= static_cast<entry_t>(-1.5);
  cout << "Matrix M = " << endl << M << endl << endl;
}
/* SAM_LISTING_END_4 */

int main() {
  MatrixXd M(6, 7);
  // Fill matrix by accessing entries directly
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < M.cols(); j++) {
      M(i, j) = i - j;
    } 
  }

  cout << "blockAccess(M)" << endl << endl;
  blockAccess(M);
}