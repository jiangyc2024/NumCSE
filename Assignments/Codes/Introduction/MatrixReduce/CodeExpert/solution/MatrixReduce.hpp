#ifndef MATRIXREDUCE_HPP
#define MATRIXREDUCE_HPP

// The purpose of this exercise is introduce reduction operations.

// TO DO: Include the appropriate header files and namespaces
// START
#include <Eigen/Dense>
#include <iostream>
using namespace Eigen;
// END

// Eigen matrices have a range of methods that reduce them to a single number.
// Suppose A is a MatrixXd object, then a few useful methods are:
// A.size(): the number of entries in A,
// A.sum(): the sum of all entries in A,
// A.prod(): the product of all entries in A,
// A.minCoeff(): the minimum value amongst the entries of A,
// A.maxCoeff(): the maximum value amongst the entries of A, and
// A.norm(): The (Frobenius) norm of the matrix A.

// TO DO: Write a function "average" that takes as argument a "MatrixXd" object
// and returns a "double". m = average(A) should calculate the average of the
// entries of A.
// START
/* SAM_LISTING_BEGIN_0 */
double average(const MatrixXd &A) {
  /*
   * Calculates the average of the entries of A
   */
  double m;
  m = A.sum() / A.size();
  return m;
}
/* SAM_LISTING_END_0 */
// END

// Question: is there a reduction method that returns the average of A directly?

/* SAM_LISTING_BEGIN_1 */
double percent_zero(const MatrixXd &A) {
  /*
   * Calculates how many entries of A are exactly equal to zero,
   * as a percentage of the total number of entries.
   */

  // Eigen provides an Array class, that is similar to Matrix,
  // but has operations that are entry-wise rather than linear algebra
  // operations. For example, multiplying two Arrays is entry-wise
  // multiplication, and not matrix multiplication.
  // The benefit of using an Array here is that we can do entry-wise boolean
  // comparison, which is not available for matrices.
  ArrayXXd Arr = A.array();

  double ratio;
  // TO DO: Calculate the ratio of zero entries in A.
  // START
  int zeros = (Arr == 0).count();
  ratio = ((double)zeros) / A.size();
  // END

  return ratio * 100;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
bool has_zero_column(const MatrixXd &A) {
  /*
   * Returns 1 if A has a column that is exactly equal to zero.
   * Returns 0 otherwise.
   */

  bool result;

  // A vector is the zero vector if and only if it has norm 0.
  // The following vector contains the squared norms of the columns of A.
  VectorXd norms = A.colwise().squaredNorm();
  // The norm is 0 if and only if the norm squared is zero.
  // We use squaredNorm() instead of norm(), because norm() =
  // sqrt(squaredNorm()) calculates a square root that is not necessary for our
  // purposes.

  // TO DO: Check if any one of the norms is equal to zero.
  // Hint: Use an array to perform entry-wise comparison.
  // START
  result = (norms.array() == 0).any();
  // END

  return result;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
MatrixXd columns_sum_to_zero(const MatrixXd &A) {
  /*
   * Returns a matrix that is like A, but the entries on the diagonal
   * have been changed so that the sum of the columns is equal to 0.
   */

  MatrixXd B(A);

  // TO DO: Replace the diagonal of B with values such that the columns of B sum
  // up to zero. Hint: Use diagonal(), rowwise(), and sum().
  // START
  int p = std::min(B.rows(), B.cols());
  B.diagonal() = VectorXd::Zero(p);
  B.diagonal() = -B.rowwise().sum();
  // END

  return B;
}

/* SAM_LISTING_END_3 */

#endif
