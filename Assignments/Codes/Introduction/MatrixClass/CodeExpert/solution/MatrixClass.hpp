#ifndef MATRIXCLASS_HPP
#define MATRIXCLASS_HPP

// The purpose of this exercise is to introduce the Eigen::Matrix class.

#include <iostream>

// TO DO: Include the "Eigen/Dense" header file with #include <Eigen/Dense>
// START
#include <Eigen/Dense>
// END

/* SAM_LISTING_BEGIN_0 */
Eigen::Matrix<double, 2, 2> smallTriangular(double a, double b, double c) {
  /*
   * This functions returns a 2 by 2 triangular matrix of doubles a, b, c.
   */

  // We know in advance that we want to create a matrix of doubles with a fixed
  // size of 2 by 2. Therefore, we pass the parameters <double, 2, 2> to the
  // template class Eigen::Matrix.
  Eigen::Matrix<double, 2, 2> A;

  // We have declared the variable A of the type Eigen::Matrix<double,2,2>,
  // but we have not initialized its entries.
  // We can do this using comma-initialization:
  A << a, b, 0, c;

  // Question: Is A now an upper triangular matrix, or a lower triangular
  // matrix?

  return A;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd constantTriangular(int n, double val) {
  /*
   * This function returns an n by n upper triangular matrix with the constant
   * value val in the upper triangular part.
   */

  // Now we do not know the size of our matrix at compile time.
  // Hence, we use the special value Eigen::Dynamic to set the size of A.
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A;

  // Eigen::Matrix has a method Zero(n_rows,n_cols) which returns the n_rows by
  // n_cols matrix whose entries are all equal to zero.
  A = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);

  // To get a triangular matrix, we must set the values above the diagonal to
  // val. This we can do by a nested for-loop.
  for (int row = 0; row < n; row++) {

    // TO DO: Write the inner for-loop.
    // Hint: We can access and change the entries of A using parentheses, i.e.
    // A(row,col) = val; Note that indexing starts at the value 0 (as usual),
    // and not 1 (as in Matlab).

    // START
    for (int col = row; col < n; col++) {
      A(row, col) = val;
    }
    // END
  }
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double casting() {
  /*
   * This function does not do anything meaningful.
   * It is only meant to introduce vectors and how to typecast.
   */

  // Because the syntax Eigen::Matrix< type, n_rows, n_cols > is very
  // cumbersome, Eigen provides convenience classes that work as shorthands.
  // For example, Eigen::Matrix2d is shorthand for Eigen::Matrix< double, 2, 2
  // >.

  // Vectors are just a special type of matrices that have only one column.
  // Thus, VectorXi is a shorthand for Eigen::Matrix< int, Eigen::Dynamic, 1 >.
  // Constant(2,1) creates a vector of size 2 and initializes the entries with
  // the value 1.
  Eigen::VectorXi u = Eigen::VectorXi::Constant(2, 1);

  // std::complex is a template class for complex numbers.
  // Here we declare two complex numbers, with real and imaginary parts
  // represented by doubles. z0 = 1 - i z1 = 5 + i
  std::complex<double> z0(1, -1), z1(5, 1);

  // Declare and initialize a size 2 vector of std::complex<double>.
  Eigen::Vector2cd v;
  v(0) = z0;
  v(1) = z1;

  double x = 0.0;
  // TO DO: Calculate the inner product of u and v, and store the result in x.
  // Hint: The inner product of two vectors is given by u.dot(v), but
  // first we need to cast the "int" vector u to a "std::complex<double>"
  // vector. Use u.cast< NEW TYPE >() to achieve this. The result of the inner
  // product will be 1*(1-i) + 1*(5+i) = 6 + 0i, a real number. You can get the
  // real part of an std::complex<double> using the method "real()".
  // START
  std::complex<double> z = u.cast<std::complex<double>>().dot(v);
  x = z.real();
  // END
  
  return x;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXcd arithmetics(int n) {
  /*
   * This function does not do anything meaningful.
   * It is only meant to show some simple Eigen arithmetics and matrix
   * initializations.
   */

  // This declares dynamic size (signified by the letter 'X') matrices of
  // complex numbers (signified by the letter c) with real and imaginary
  // parts represented by doubles (signified by 'd').
  Eigen::MatrixXcd A, C, I;

  // We initialize the matrices arbitrarily:
  // an n by n lower triangular matrix,
  A = constantTriangular(n, 5.0).transpose().cast<std::complex<double>>();
  // The n by n identity matrix,
  I = Eigen::MatrixXcd::Identity(n, n);

  // Declare the n by n matrix B.
  Eigen::MatrixXcd B(n, n);

  // TO DO: Fill in the matrix B such that B(k,l) = (k+l*i)/(k-l*i),
  // where i is the imaginary unity, i^2 = -1.
  // START
  for (int k = 0; k < n; k++) {
    for (int l = 0; l < n; l++) {
      std::complex<double> tmp(k + 1, l + 1);
      B(k, l) = tmp / std::conj(tmp);
    }
  }
  // END

  // We can perform arithmetics on matrices: +, -, *
  // Note that for + and -, the left hand matrix and the right hand matrix have
  // to be the same size, and that matrix multiplication * is only defined if
  // the number of columns in the left hand matrix is the same as the number of
  // rows in the right hand matrix.
  C = B * (A - 5.0 * I);

  // Eigen::VectorXcd is shorthand for
  // Eigen::Matrix< std::complex<double>, Eigen::Dynamic, 1>.
  // Hence, all the same arithmetics work for vectors.
  Eigen::VectorXcd u, v;

  // TO DO: Initialize the entries of u as the integers from 1 to n,
  // multiply u by the matrix C, and store the result in v.
  // START
  u = Eigen::VectorXcd::LinSpaced(n, 1, n); // Or use a for-loop.
  v = C * u;
  // END

  return v;
}
/* SAM_LISTING_END_3 */

#endif
