#include <iostream>
#include <iomanip>
#include <math.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "timer.h"

using namespace Eigen;

/* \brief Compute powers of a square matrix
 * Use smart binary representation
 * \param[in,out] A square (complex) matrix. $A^k$ is stored in $A$
 * \param[in] k positive integer for $A^k$
 */
/* SAM_LISTING_BEGIN_1 */
MatrixXcd matPow(MatrixXcd & A, unsigned int k) {
    // Identity matrix
    MatrixXcd X = MatrixXcd::Identity(A.rows(), A.cols());

    // $p$ is used as binary mask to check wether
    // given $k = \sum_{i = 0}^M b_i 2^i$, $b_i = 1$
    // (i.e. if $k$ has a 1 in the $i$-th binary digit)
    // obtaining the binay representation of p can be done in many ways,
    // here we use $\tilde{}k \& p$ to check if $b_i = 1$
    unsigned int p = 1;

    // Cycle all the way up to the 1st 1 in the binary
    // representation of $k$
    for(unsigned int j = 1; j <= ceil(log2(k)); ++j) {

      if( ( ~k & p ) == 0 ) { // if( $b_i != 0$ )
            X = X*A;
        }


        A = A*A;
        p = p << 1; // p = p*2;
    }
    return X;
}
/* SAM_LISTING_END_1 */

/* \brief Creates a Vandermonde matrix with size $n$
 * \param[in] n size of the output matrix
 * \return Vandermonde matrix with entries: $exp(2\pi i j k / n) / sqrt(n)$
 */
MatrixXcd construct_matrix(unsigned int n) {

    // This is $\pi$
    double PI = M_PI; // from math.h
    // This is $i$, the imaginary unit
    std::complex<double> I = std::complex<double>(0,1);

    // Allocate space for $\mathbf{A}$
    MatrixXcd A(n,n);

    // Fill matrix
    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
          A(i,j) = std::exp(2. * PI * I *
                            (double)(i+1) * (double)(j+1) / (double) n)
                 / std::sqrt((double) n);
        }
    }

    return A;
}

int main(void) {
    // Check/Test with provided, complex, matrix
    unsigned int n = 3; // size of matrix
    unsigned int k = 9; // power

    MatrixXcd A = construct_matrix(n);

    // Output results
    MatrixXcd Eigen_A_pow_k = A.pow(k);
    std::cout << "A = " << std::endl
              << A << std::endl
              << "-> Eigen:" << std::endl
              << "A^" << k << " = " << std::endl
              << Eigen_A_pow_k << std::endl;
    MatrixXcd A_pow_k = matPow(A, k);
    std::cout << "-> Ours:"  << std::endl
              << "A^" << k << " = " << std::endl
              << A_pow_k << std::endl;

    // Use this to test your implementation, error should be very small (machine precision)
    std::cout << "-- Error of our implementation VS Eigen: "
              << (A_pow_k - Eigen_A_pow_k).norm() << std::endl;

    /* SAM_LISTING_BEGIN_2 */

    // Will contain runtimes for each k, k stored in powers
    std::vector<double> powers, times_matPow, times_Eigen_pow;

    std::cout << "-- Table of runtimes" << std::endl;
    std::cout << std::setw(15) << "K"
              << std::setw(15) << "Own matPow"
              << std::setw(15) << "Eigen pow()" << std::endl;
    // Loop from $2$ to $2^31$
    for(unsigned int K = 2; K <= (1u << 31); K = K << 1) {

        unsigned int repeats = 10;
        Timer tm_pow, tm_Eigen_pow;

        // Repeat the test repeat times
        for(unsigned int r = 0; r < repeats; ++r) {

            // Build Vandermonde matrix of size n
            MatrixXcd X, A = construct_matrix(n);

            // Compute runtime if own implementation
            tm_pow.start();
            X = matPow(A, K);
            tm_pow.stop();

            // Compute runtime of eigen implementation
            A = construct_matrix(n);
            tm_Eigen_pow.start();
            X = A.pow(K);
            tm_Eigen_pow.stop();
        }

        // Output table
        std::cout << std::setw(15) << K
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << tm_pow.min()
                  << std::setw(15) << tm_Eigen_pow.min() << std::endl;

        // Store data
        powers.push_back( K );
        times_matPow.push_back( tm_pow.min() );
        times_Eigen_pow.push_back( tm_Eigen_pow.min() );
    }


    /* SAM_LISTING_END_2 */
}
