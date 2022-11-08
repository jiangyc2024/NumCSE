#ifndef MATRIXBLOCK_HPP
#define MATRIXBLOCK_HPP

// The purpose of this exercise is introduce block operations on Eigen matrices.

#include <cassert>

// TODO: Include the appropriate header files.
// START
#include <Eigen/Dense>
// END

/**
 * @brief This function takes a matrix A, and returns a matrix that is exactly
 * the same, except that row p and column q are set to zero.
 *
 * @param A The matrix
 * @param p Row to set to zero
 * @param q Column to set to zero
 * @return Modified A
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd zero_row_col(const Eigen::MatrixXd& A, unsigned int p,
                             unsigned int q) {
  // assert that we do not run into Eigen assertions.
  assert(A.rows() > p && A.cols() > q && "p or q out of bounds.");

  // Make a copy of A.
  Eigen::MatrixXd Anew(A);

  // TODO: (0-2.a) Set the entries of row number p and column number q to zero.
  // Hint: We can access rows and columns of Anew by Anew.row() and Anew.col().
  // The method setZero() is useful here.
  // START
  Anew.row(p).setZero();
  Anew.col(q).setZero();
  // END

  return Anew;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Writing as a block matrix A = [B C], where B denotes the first p
 *        columns of A, and C denotes the q=(A.cols() - p) last columns, this
 *        functions returns D = [C B].
 *
 * @param A The matrix whose columns should be swapped
 * @param p The amount of columns to swap
 * @return The swapped matrix
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd swap_left_right_blocks(const Eigen::MatrixXd& A,
                                       unsigned int p) {
  assert(A.cols() >= p && "Can only swap at most A.cols() columns.");
  Eigen::MatrixXd B, C;

  // We can use .rows() and .cols() to get the number of rows and columns in A.
  unsigned int q = A.cols() - p;

  // Make a copy of A.
  Eigen::MatrixXd Anew(A);

  // TODO: (0-2.b) Swap the first p columns of A with the rightmost q columns
  // and write them in Anew. The methods block(), leftCols()/rightCols() may
  // come in handy.
  // Hint: A.block( i, j, m, n ) returns the m by n block that has its top-left
  // corner at the index (i,j) in A.
  // START
  Anew.leftCols(q) = A.block(0, p, A.rows(), q);
  Anew.rightCols(p) = A.block(0, 0, A.rows(), p);
  // END

  // Tip: Many more methods exist that are special cases of block(),
  // e.g. topRows(), bottomRows(), topLeftCorner(), bottomLeftCorner(), ...
  // For vectors we have head(), tail(), and segment().

  return Anew;
}
/* SAM_LISTING_END_1 */

/**
 * @brief This function creates an n by n tridiagonal matrix with the values
 *        a on the first subdiagonal,
 *        b on the diagonal, and
 *        c on the first superdiagonal.
 *
 *        Example for n=5:
 *            [ b c 0 0 0 ]
 *            [ a b c 0 0 ]
 *        A = [ 0 a b c 0 ]
 *            [ 0 0 a b c ]
 *            [ 0 0 0 a b ]
 *
 * @param n Dimension
 * @param a Value of first subdiagonal
 * @param b Value of diagonal
 * @param c Value of first superdiagonal
 * @return Tridiagonal matrix
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::MatrixXd tridiagonal(unsigned int n, double a, double b, double c) {
  Eigen::MatrixXd A;

  // TODO: (0-2.c) Fill the matrix A with zeros. Then fill the subdiagonal with
  // the value a, the diagonal with b and the superdiagonal with c.
  // Hint: You can get the diagonal of A by A.diagonal(). Moreover, you can get
  // the super- and subdiagonals by passing +1 or -1 as arguments to
  // A.diagonal(). Here, however, we can use templated access because the
  // distance to the diagonal is known at compile time. This method is usually
  // faster.
  // START
  A = Eigen::MatrixXd::Zero(n, n);
  A.diagonal<-1>() = Eigen::VectorXd::Constant(n - 1, a);
  A.diagonal() = Eigen::VectorXd::Constant(n, b);
  A.diagonal<1>() = Eigen::VectorXd::Constant(n - 1, c);
  // END

  return A;
}
/* SAM_LISTING_END_2 */

#endif
