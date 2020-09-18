////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include "timer.h"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace Eigen;

/* \brief Compute the Matrix product $A \times B$ using Strassen's algorithm.
 * \param[in] A Matrix $2^k \times 2^k$
 * \param[in] B Matrix $2^k \times 2^k$
 * \param[out] Matrix product of A and B of dim $2^k \times 2^k$
 */
/* SAM_LISTING_BEGIN_0 */
MatrixXd strassenMatMult(const MatrixXd &A, const MatrixXd &B) {
  // Ensure square matrix
  assert(A.rows() == A.cols() && "Matrix A must be square");
  assert(B.rows() == B.cols() && "Matrix B must be square");
  // Matrix dimension must be a power of 2
  assert(A.rows() % 2 == 0 && "Matrix dimensions must be a power of two.");

  const unsigned n = A.rows();

  // The function is recursive and acts on blocks of size $n/2 \times n/2$
  // i.e. exploits fast product of 2x2 block matrix
  if (n == 2) { // End of recursion
    MatrixXd C(2, 2);
    C << A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0),
        A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1),
        A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0),
        A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);
    return C;
  }

  MatrixXd Q0(n / 2, n / 2), Q1(n / 2, n / 2), Q2(n / 2, n / 2),
      Q3(n / 2, n / 2), Q4(n / 2, n / 2), Q5(n / 2, n / 2), Q6(n / 2, n / 2);

  MatrixXd C(n, n);
  // TO DO: (2-4.a) Finish Strassen's algorithm.
  // Hint: Use strassenMatMult() to fill in the matrices Q0,..., Q6,
  // and then use Q0,..., Q6 to fill in the matrix C.
  // Note that comma-initialization works for block matrices.
  // START:
  MatrixXd A11 = A.topLeftCorner(n / 2, n / 2);
  MatrixXd A12 = A.topRightCorner(n / 2, n / 2);
  MatrixXd A21 = A.bottomLeftCorner(n / 2, n / 2);
  MatrixXd A22 = A.bottomRightCorner(n / 2, n / 2);

  MatrixXd B11 = B.topLeftCorner(n / 2, n / 2);
  MatrixXd B12 = B.topRightCorner(n / 2, n / 2);
  MatrixXd B21 = B.bottomLeftCorner(n / 2, n / 2);
  MatrixXd B22 = B.bottomRightCorner(n / 2, n / 2);

  Q0 = strassenMatMult(A11 + A22, B11 + B22);
  Q1 = strassenMatMult(A21 + A22, B11);
  Q2 = strassenMatMult(A11, B12 - B22);
  Q3 = strassenMatMult(A22, B21 - B11);
  Q4 = strassenMatMult(A11 + A12, B22);
  Q5 = strassenMatMult(A21 - A11, B11 + B12);
  Q6 = strassenMatMult(A12 - A22, B21 + B22);

  C << Q0 + Q3 - Q4 + Q6, Q2 + Q4, Q1 + Q3, Q0 + Q2 - Q1 + Q5;
  // END

  return C;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
double test_strassen() {
  // Check algorithm for correctness

  // Size of the matrix is a power of 2.
  int k = 4;
  int n = std::pow(2, k);
  double err = 1;

  // TO DO: (2-4.b) Generate two random matrices using MatrixXd::Random(n,n).
  // and compare the result of strassenMatMult() and the built-in matrix
  // multiplication. Save the difference in the double err. START
  MatrixXd A = MatrixXd::Random(n, n);
  MatrixXd B = MatrixXd::Random(n, n);

  // Testing matrix multiplication
  MatrixXd AB = strassenMatMult(A, B);
  // Eigen Matrix multiplication
  MatrixXd AxB = A * B;

  err = (AB - AxB).norm();
  // END

  return err;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void time_strassen() {
  // Unfortunately, the time limit set by Code Expert may not be sufficiently
  // long to show the advantage of Strassen's algorithm for large matrices.

  // Minimum number of repetitions
  unsigned int repetitions = 5;

  // Display header column
  std::cout << std::setw(4) << "k" << std::setw(15) << "A*B" << std::setw(15)
            << "Strassen" << std::endl;
  for (unsigned k = 3; k <= 6; k++) {
    unsigned int n = std::pow(2, k);

    // Initialize random input matricies
    MatrixXd A = MatrixXd::Random(n, n);
    MatrixXd B = MatrixXd::Random(n, n);

    // Timer to collect runtime of each individual run
    Timer timer, timer_own;

    // TO DO: Use the .start(), and .stop() methods of timer and timer\_own to
    // measure the runtimes of Eigen's built-in multiplication and Strassen's
    // algorithm. Perform a few trials for each k (use the variable
    // repetitions). START
    for (unsigned int r = 0; r < repetitions; ++r) {
      // initialize memory for result matrix
      MatrixXd AxB, AxB2;

      // Benchmark Eigen's matrix multiplication
      timer.start(); // start timer
      AxB = (A * B); // do the multiplication
      timer.stop();  // stop timer

      // Benchmark Strassen's matrix multiplication
      timer_own.start();            // start timer
      AxB2 = strassenMatMult(A, B); // do the multiplication
      timer_own.stop();             // stop timer
    }
    // END

    // Print runtimes
    std::cout << std::setw(4) << k // power
              << std::setprecision(3) << std::setw(15) << std::scientific
              << timer.min() // eigen timing
              << std::setprecision(3) << std::setw(15) << std::scientific
              << timer_own.min() // strassing timing
              << std::endl;
  }
}
/* SAM_LISTING_END_2 */
