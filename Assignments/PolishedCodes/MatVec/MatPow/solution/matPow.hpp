#ifndef MATPOW_HPP
#define MATPOW_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>

#include "timer.h"

/**
 * \brief Compute powers of a square matrix. Use a smart binary representation.
 *
 * \param A square (complex) matrix.
 * \param k Positive integer for $A^k$
 * \return Eigen::MatrixXcd the matrix exponential
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXcd matPow(Eigen::MatrixXcd& A, unsigned int k) {
  // Identity matrix
  Eigen::MatrixXcd X = Eigen::MatrixXcd::Identity(A.rows(), A.cols());

  // TODO: (1-6.a) Compute the k-th power of A in X.
  // START

  // $p$ is used as binary mask to check whether,
  // given $k = \sum_{i = 0}^M b_i 2^i$, $b_i = 1$
  // (i.e. if $k$ has a 1 in the $i$-th binary digit).
  // Obtaining the binary representation of p can be done in many ways,
  // here we use $\sim k \& p$ to check if $b_i = 1$
  unsigned int p = 1;

  // Cycle all the way up to the 1st 1 in the binary
  // representation of $k$ (this is called the "most significant bit")
  for (unsigned int j = 1; j <= std::ceil(std::log2(k)); ++j) {
    if ((~k & p) == 0) {  // if( $b_i \neq 0$ )
      X = X * A;
    }
    A = A * A;
    p = p << 1;  // p = p*2;
  }
  // END
  return X;
}
/* SAM_LISTING_END_1 */

/**
 * \brief Creates a Vandermonde matrix with size $n$
 *
 * \param n size of the output matrix
 * \return Eigen::MatrixXcd Vandermonde matrix with entries: $exp(2\pi i j k /
 * n) / sqrt(n)$
 */
Eigen::MatrixXcd construct_matrix(unsigned int n) {
  // This is $\pi$
  constexpr double PI = M_PI;  // from math.h
  // This is $i$, the imaginary unit
  std::complex<double> I = std::complex<double>(0, 1);

  // Allocate space for $\mathbf{A}$
  Eigen::MatrixXcd A(n, n);

  // Fill matrix
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      A(i, j) = std::exp(2. * PI * I * (double)(i + 1) * (double)(j + 1) /
                         (double)n) /
                std::sqrt((double)n);
    }
  }
  return A;
}

/* SAM_LISTING_BEGIN_2 */
void tabulateRuntime(unsigned int n) {
  std::cout << "-- Table of runtimes" << std::endl;
  std::cout << std::setw(15) << "K" << std::setw(15) << "Own matPow"
            << std::setw(15) << "Eigen pow()" << std::endl;
  constexpr unsigned int repeats = 10;
  // Loop from $2$ to $2^31$
  for (unsigned long long int K = 2; K <= (1u << 31); K = K << 1) {
    Timer tm_pow, tm_Eigen_pow;

    // Repeat the test repeat times
    for (unsigned int r = 0; r < repeats; ++r) {
      // TODO: (1-6.c) Measure runtime of your own and Eigen's implementation.
      // You may use construct_matrix.
      // START
      // Build Vandermonde matrix of size n
      Eigen::MatrixXcd X, A = construct_matrix(n);

      // Compute runtime if own implementation
      tm_pow.start();
      X = matPow(A, K);
      tm_pow.stop();

      // Compute runtime of eigen implementation
      A = construct_matrix(n);
      tm_Eigen_pow.start();
      X = A.pow(K);
      tm_Eigen_pow.stop();
      // END
    }

    // Output table
    std::cout << std::setw(15) << K << std::scientific << std::setprecision(3)
              << std::setw(15) << tm_pow.min() << std::setw(15)
              << tm_Eigen_pow.min() << std::endl;
  }
}
/* SAM_LISTING_END_2 */

#endif