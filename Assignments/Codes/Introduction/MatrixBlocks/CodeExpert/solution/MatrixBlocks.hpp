#ifndef MATRIXBLOCK_HPP
#define MATRIXBLOCK_HPP

// The purpose of this exercise is introduce block operations on Eigen matrices.

// TO DO: Include the appropriate header files.
// START
#include <Eigen/Dense>
// END

// We can use the Eigen namespace to improve the readability of our code.
// This allows us to skip the "Eigen::" in "Eigen::MatrixXd" for example.
// TO DO: Add the following line below: using namespace Eigen;
// START
using namespace Eigen;
// END

/* SAM_LISTING_BEGIN_0 */
MatrixXd zero_row_col(const MatrixXd &A, int p, int q) {
  /*
   * This function takes a matrix A, and returns a matrix that is exactly the
   * same, except that row p and column q are set to zero.
   */
  // Make a copy of A.
  MatrixXd Anew(A);

  // TO DO: Set the entries of row number p and column number q to zero.
  // Hint: We can access rows and columns of A by A.row() and A.col().
  // The setZero() is useful here. START
  Anew.row(p).setZero();
  Anew.col(q).setZero();
  // END

  return Anew;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
MatrixXd swap_left_right_blocks(const MatrixXd &A, int p) {
  /*
   * Writing as a block matrix A = [B C], where B denotes the first p columns of
   * A, and C denotes the q=(A.cols() - p) last columns, this functions returns
   * D = [C B].
   */

  MatrixXd B, C;

  // We can use .rows() and .cols() to get the number of rows and columns in A.
  int q = A.cols() - p;

  // A.block( i, j, m, n ) returns the m by n block that has its top-left corner
  // at the index (i,j) in A. Hence, the first p columns of A can be accessed in
  // the following way:
  B = A.block(0, 0, A.rows(), p);

  // TO DO: Use A.block() to define C as the matrix containing the last q
  // columns of A. START
  C = A.block(0, p, A.rows(), q);
  // END

  // Make a copy of A.
  MatrixXd Anew(A);
  // The block() method can access arbitrary blocks within a matrix.
  // For our purposes, it is actually simpler to use leftCols() and rightCols().
  Anew.leftCols(q) = C;

  // TO DO: Use A.rightCols() to fill in the remaining columns of the new matrix
  // A. START
  Anew.rightCols(p) = B;
  // END

  // Tip: Many more methods exist that are special cases of block(),
  // e.g. topRows(), bottomRows(), topLeftCorner(), bottomLeftCorner(), ...
  // For vectors we have head(), tail(), and segment().

  return Anew;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
MatrixXd tridiagonal(int n, double a, double b, double c) {
  /*
   * This function creates an n by n tridiagonal matrix with the values
   *      a on the first subdiagonal,
   *      b on the diagonal, and
   *      c on the first superdiagonal.
   * Example for n=5:
   *     [ b c 0 0 0 ]
   *     [ a b c 0 0 ]
   * A = [ 0 a b c 0 ]
   *     [ 0 0 a b c ]
   *     [ 0 0 0 a b ]
   */

  MatrixXd A;

  // TO DO: Fill the matrix A with zeros. Then fill the subdiagonal with the
  // value a, the diagonal with b and the superdiagonal with c. Hint: You can
  // get the diagonal of A by A.diagonal(). Moreover, you can get the super- and
  // subdiagonals by passing +1 or -1 as arguments to A.diagonal(). START
  A = MatrixXd::Zero(n, n);
  A.diagonal(-1) = VectorXd::Constant(n - 1, a);
  A.diagonal() = VectorXd::Constant(n, b);
  A.diagonal(1) = VectorXd::Constant(n - 1, c);
  // END
  return A;
}
/* SAM_LISTING_END_2 */

#endif
