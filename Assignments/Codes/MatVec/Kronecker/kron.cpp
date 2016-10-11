#include <iostream>
#include <iomanip>

#include <vector>

#include <Eigen/Dense>

#if INTERNAL
#include <figure/figure.hpp>
#endif // INTERNAL
#include "timer.h"

using namespace Eigen;

/* \brief Compute the Kronecker product.
 * Computes $\mathbf{C} = \mathbf{A} \otimes \mathbf{B}$.
 * \param[in] A Matrix of size $n \times n$
 * \param[in] B Matrix of size $n \times n$
 * \param[out] C Kronecker product of A and B of dim $n^2 \times n^2$
 */
/* SAM_LISTING_BEGIN_1 */
void kron(const MatrixXd & A, const MatrixXd & B,
          MatrixXd & C) {
#if SOLUTION
    // Allocate enough space for the matrix
    C = MatrixXd(A.rows()*B.rows(), A.cols()*B.cols());
    for(unsigned int i = 0; i < A.rows(); ++i) {
        for(unsigned int j = 0; j < A.cols(); ++j) {
            // We use eigen block operations to set the values of
            // each $n \times n$ block.
            C.block(i*B.rows(),j*B.cols(), B.rows(), B.cols())
                = A(i,j)*B; // $\in \mathbb{R}^{(n \times n)}$
        }
    }
#else // TEMPLATE
    // TODO: Implement the Kronecker product
    // Hint: Use MatrixXd::block(int, int, int, int)
    // See https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
#endif // TEMPLATE
}
/* SAM_LISTING_END_1 */

/* \brief Compute the Kronecker product applied to a vector.
 * Computes $\mathbf{y} = (\mathbf{A} \otimes \mathbf{B}) \mathbf{x}$.
 * \param[in] A Matrix of size $n \times n$
 * \param[in] B Matrix of size $n \times n$
 * \param[in] x Vector of dim. $n^2 \times n^2$
 * \param[out] y Vector y = kron(A,B)*x
 */
/* SAM_LISTING_BEGIN_2 */
void kron_mult(const MatrixXd &A, const MatrixXd &B,
               const VectorXd &x, VectorXd &y) {
    assert(A.rows() == A.cols() &&
           A.rows() == B.rows() &&
           B.rows() == B.cols() &&
           "Matrices A and B must be square matrices with same size!");
    assert(x.size() == A.cols()*A.cols() &&
           "Vector x must have length A.cols()^2");
    unsigned int n = A.rows();

#if SOLUTION
    // Allocate space for output
    y = VectorXd::Zero(n*n);

    // Note: this is like doing a matrix-vector multiplication
    // where the entries of the matrix are smaller matrices
    // and entries of the vector are smaller vectors

    // Loop over all segments of x ($\tilde{x}$)
    for(unsigned int j = 0; j < n; ++j) {
        // Reuse computation of z
        VectorXd z = B * x.segment(j*n, n);
        // Loop over all segments of y
        for(unsigned int i = 0; i < n; ++i) {
            y.segment(i*n, n) += A(i, j)*z;
        }
    }
#else // TEMPLATE
    // TODO: implement $y = (A \otimes B) \cdot x$ as efficiently as possible
#endif // TEMPLATE
}
/* SAM_LISTING_END_2 */

/* \brief Compute the Kronecker product $C = A \otimes B$.
 * Use fast reshaping (similar to Matlab reshape)
 * WARNING: using Matrix::Map we assume the matrix is in Column major format,
 *          the code is not valid for Row Major format.
 * \param[in] A Matrix $n \times n$
 * \param[in] B Matrix $n \times n$
 * \param[in] x Vector of dim $n^2$
 * \param[out] y Vector y = kron(A,B)*x
 */
/* SAM_LISTING_BEGIN_3 */
void kron_reshape(const MatrixXd & A, const MatrixXd & B,
                  const VectorXd & x, VectorXd & y) {
    assert(A.rows() == A.cols() && A.rows() == B.rows() && B.rows() == B.cols() &&
           "Matrices A and B must be square matrices with same size!");
    unsigned int n = A.rows();

#if SOLUTION
    MatrixXd t = B * MatrixXd::Map(x.data(), n, n) * A.transpose();
    y = MatrixXd::Map(t.data(), n*n, 1);
#else // TEMPLATE
    // TODO: compute y
    // Hint: Use MatrixXd::Map(...)
    // See https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
#endif
}
/* SAM_LISTING_END_3 */

int main(void) {

    // Testing correctness of Kron
    MatrixXd A(2,2);
    A << 1, 2, 3, 4;
    MatrixXd B(2,2);
    B << 5, 6, 7, 8;
    MatrixXd C;

    VectorXd x = Eigen::VectorXd::Random(4);
    VectorXd y;

    std::cout << "Testing kron, kron_mult, and kron_reshape with small matrices..."
              << std::endl;

    // Compute using kron
    kron(A, B, C);
    y = C*x;

    std::cout << "kron(A,B) = " << std::endl
              << C << std::endl;
    std::cout << "Using kron: y = " << std::endl
              << y << std::endl;

    // Compute using kron_mult
    kron_mult(A, B, x, y);
    std::cout << "Using kron_mult: y =" << std::endl
              << y << std::endl;

    // Compute using kron_reshape
    kron_reshape(A, B, x, y);
    std::cout << "Using kron_reshape: y =" << std::endl
              << y << std::endl;

#if SOLUTION
    // Compute runtime of different implementations of kron

    /* SAM_LISTING_BEGIN_4 */
    // We repeat each runtime measurment 10 times
    // (this is done in order to remove outliers)
    unsigned int repeats = 10;
#if INTERNAL
    std::vector<double> sizes, times_kron, times_kron_mult, times_kron_map;
#endif // INTERNAL

    std::cout << "Runtime for each implementation." << std::endl;
    std::cout << std::setw(5) << "n"
              << std::setw(15) << "kron"
              << std::setw(15) << "kron_mult"
              << std::setw(15) << "kron_reshape"
              << std::endl;
    // Loop from $M = 2,\dots,2^8$
    for(unsigned int M = 2; M <= (1 << 8); M = M << 1) {
        Timer tm_kron, tm_kron_mult, tm_kron_map;
        // Run experiments "repeats" times
        for(unsigned int r = 0; r < repeats; ++r) {
            // Random matrices for testing
            A = MatrixXd::Random(M,M);
            B = MatrixXd::Random(M,M);
            x = VectorXd::Random(M*M);

            // Do not want to use kron for large values of M
            if( M < (1 << 6)) {
                // Kron using direct implementation
                tm_kron.start();
                kron(A,B,C);
                y = C*x;
                tm_kron.stop();
            }

            // Kron matrix-vector multiplication
            tm_kron_mult.start();
            kron_mult(A,B,x,y);
            tm_kron_mult.stop();

            // Kron using reshape
            tm_kron_map.start();
            kron_reshape(A,B,x,y);
            tm_kron_map.stop();
        }

        double kron_time = (M < (1 << 6)) ? tm_kron.min() : std::nan("");
        std::cout << std::setw(5) << M
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << kron_time
                  << std::setw(15) << tm_kron_mult.min()
                  << std::setw(15) << tm_kron_map.min() << std::endl;
#if INTERNAL
        sizes.push_back( M );
        times_kron.push_back( kron_time );
        times_kron_mult.push_back( tm_kron_mult.min() );
        times_kron_map.push_back( tm_kron_map.min() );
#endif // INTERNAL

    }
    /* SAM_LISTING_END_4 */
#else // TEMPLATE
    // TODO: Compute runtime of kron, kron_mult and kron_reshape.
    // Output the runtimes in scientific notation with 3 decimal digits
    // Use the Timer class.
#endif // TEMPLATE

#if INTERNAL
  mgl::Figure fig;
  fig.title("Timings of kron");
  fig.ranges(2, 9000, 1e-8, 1e3);
  fig.setlog(true, true); // set loglog scale
  fig.plot(sizes, times_kron, " r+").label("runtime");
  fig.plot(sizes, times_kron_mult, " b+").label("runtime");
  fig.plot(sizes, times_kron_map, " g+").label("runtime");
  fig.fplot("1e-9*x^4", "k|").label("O(n^4)");
  fig.xlabel("Matrix size (M))");
  fig.ylabel("Time [s]");
  fig.legend(0, 1);
  fig.save("kron_timing.eps");
  fig.save("kron_timing.png");
#endif // INTERNAL
}
