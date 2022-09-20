#ifndef MATRIXREDUCE_HPP
#define MATRIXREDUCE_HPP

// The purpose of this exercise is introduce reduction operations.

// TODO: Include the appropriate header files
// START

// END

// Eigen matrices have a range of methods that reduce them to a single number.
// Suppose A is a MatrixXd object, then a few useful methods are:
// A.size(): the number of entries in A,
// A.sum(): the sum of all entries in A,
// A.prod(): the product of all entries in A,
// A.minCoeff(): the minimum value amongst the entries of A,
// A.maxCoeff(): the maximum value amongst the entries of A, and
// A.norm(): The (Frobenius) norm of the matrix A.

/**
 * @brief Calculates the average of the entries of A.
 *
 * @param A matrix
 * @return double average
 */
/* SAM_LISTING_BEGIN_0 */
double average(const Eigen::MatrixXd& A) {
  double m = 0;

  // TODO: (1-3.a) Calculate the average of the entries of A.
  // START

  // END

  return m;
}
/* SAM_LISTING_END_0 */

// Question: is there a reduction method that returns the average of A directly?

/**
 * @brief Calculates how many entries of A are exactly equal to zero,
 * as a percentage of the total number of entries.
 *
 * @param A matrix
 * @return double percentage of zero entries
 */
/* SAM_LISTING_BEGIN_1 */
double percent_zero(const Eigen::MatrixXd& A) {
  // Eigen provides an Array class, that is similar to Matrix,
  // but has operations that are entry-wise rather than linear algebra
  // operations. For example, multiplying two Arrays is entry-wise
  // multiplication, and not matrix multiplication.
  // The benefit of using an Array here is that we can do entry-wise boolean
  // comparison, which is not available for matrices.
  Eigen::ArrayXXd Arr = A.array();

  double ratio = 0.;
  // TODO: (1-3.b) Calculate the ratio of zero entries in A.
  // START

  // END

  return ratio * 100;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Checks whether A has a zero column.
 *
 * @param A matrix
 * @return true if A has a column that is exactly equal to zero
 * @return false otherwise
 */
/* SAM_LISTING_BEGIN_2 */
bool has_zero_column(const Eigen::MatrixXd& A) {
  bool result;

  // A vector is the zero vector if and only if it has norm 0.
  // The following vector contains the squared norms of the columns of A.
  Eigen::VectorXd norms = A.colwise().squaredNorm();
  // The norm is 0 if and only if the norm squared is zero.
  // We use squaredNorm() instead of norm(), because norm() =
  // sqrt(squaredNorm()) calculates a square root that is not necessary for our
  // purposes.

  // TODO: (1-3.c) Check if any one of the norms is equal to zero.
  // Hint: Use an array to perform entry-wise comparison.
  // START

  // END

  return result;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Returns a matrix that is like A, but the entries on the diagonal
 * have been changed so that the sum of the columns is equal to 0.
 *
 * @param A matrix
 * @return Eigen::MatrixXd like A but with changed diagonal entries
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd columns_sum_to_zero(const Eigen::MatrixXd& A) {
  Eigen::MatrixXd B(A);

  // TODO: (1-3.d) Replace the diagonal of B with values such that the columns
  // of B sum up to zero. Hint: Use diagonal(), rowwise(), and sum().
  // START

  // END

  return B;
}
/* SAM_LISTING_END_3 */

#endif
