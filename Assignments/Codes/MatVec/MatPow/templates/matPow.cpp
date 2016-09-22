#include <iostream>
#include <iomanip>
#include <math.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "timer.h"

using namespace Eigen;

/* \brief Compute powers of a square matrix.
 * Use a smart binary representation.
 * \param[in,out] A square (complex) matrix. $A^k$ is stored in $A$
 * \param[in] k Positive integer for $A^k$
 */
/* SAM_LISTING_BEGIN_1 */
MatrixXcd matPow(MatrixXcd & A, unsigned int k) {
    // TODO: implement efficient matrix power
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

    // TODO: compute runtime of Eigen pow() and matPow().
}
