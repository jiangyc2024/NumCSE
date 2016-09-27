#include <iostream>
#include <cmath>

#include "strassen.hpp"

int main()
{
    /* SAM_LISTING_BEGIN_1 */

    // Seed the random number generator
    srand((unsigned int) time(0));

    //// Check algorithm for correctness

    // Size of the matrix
    int k = 2;
    int n = std::pow(2, k);

    // Generate random matrices
    MatrixXd A = MatrixXd::Random(n,n);
    MatrixXd B = MatrixXd::Random(n,n);

    // Testing matrix multiplication
    MatrixXd AB = strassenMatMult(A,B);
    // Eigen Matrix multiplication
    MatrixXd AxB = A*B;

    std::cout << "Using Strassen's method, A*B=" << std::endl
              << AB << std::endl;
    std::cout << "Using standard method, A*B=" << std::endl
              << AxB << std::endl;
    std::cout << "The norm of the error is: "
              << (AB-AxB).norm() << std::endl;
    /* SAM_LISTING_END_1 */
}
