#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include "timer.h"
#if INTERNAL
#include <chrono>
#include <figure/figure.hpp>
#endif // INTERNAL

using namespace Eigen;

/* @brief Compute the CCS format of matrix $A$
 * \param[in] A An $n \times n$ matrix
 * \param[out] val Vector of nonzero values of $A$
 * \param[out] row_ind Row indices of each element of 'val'
 * \param[out] col_ptr Indices of the elements in 'val' which start a column of $A$
 */
/* SAM_LISTING_BEGIN_0 */
void CCS(const MatrixXd & A, VectorXd & val, VectorXd & row_ind, VectorXd & col_ptr)
{
    // Initialization
    int n = A.size();
    MatrixXd A(n,n);

#if SOLUTION
	// Number of rows and columns
	m = A.rows();
	n = A.cols();

	// Number of nonzero entries
	int NNZ = 0;
	for(int i=0; i<m; ++i) { // Row iterator
		for(int j=0; j<n; ++j) { // Col iterator
			if(A(i,j) != 0) {
				++NNZ;
			}
		}
	}

	// Initialization
	val.resize(NNZ);
	row_ind.resize(NNZ);
	col_ptr.resize(n);
	col_ptr(0) = 0;

	// Store $A$ in CRS format
	int index = 0;
	for(int j=0; j<n; ++j) { // Col iterator
		for(int i=0; i<m; ++i) { // Row iterator
			if(A(i,j) != 0) {
				// Record the value to 'val'
				val(index) = A(i,j);
				// Record the row index to 'row_ind'
				row_ind(index) = i;
				++index;
			}
		}
		// Update 'col_ptr'
		col_ptr(j+1) = index;
	}
#else // TEMPLATE
    // TODO: compute the CCS format of matrix $A$
#endif // TEMPLATE

	return A;
}
/* SAM_LISTING_END_0 */

/* @brief Compute the CCS format of matrix $A$
 * \param[in] A An $n \times n$ matrix
 * \param[out] val Vector of nonzero values of $A$
 * \param[out] row_ind Row indices of each element of 'val'
 * \param[out] col_ptr Indices of the elements in 'val' which start a column of $A$
 */
/* SAM_LISTING_BEGIN_1 */
void CCS(const MatrixXd & A, VectorXd & val, VectorXd & row_ind, VectorXd & col_ptr)
{
    // Initialization
    int n = A.size();
    MatrixXd A(n,n);



	return A;
}
/* SAM_LISTING_END_1 */

/* @brief Solve $Ax = b$.
 * Function with naive implementation.
 * @param[in] a An $n$-dimensional vector to build the lower triangular matrix $A$
 * @param[in] b An $n$-dimensional vector
 * @param[out] x The $n$-dimensional vector $x = A^{-1}*b$
 */
/* SAM_LISTING_BEGIN_1 */
void solveA(const VectorXd & a, const VectorXd & b, VectorXd & x) {
    // Size of b, which is the size of a
    int n = b.size();
    assert( n == a.size()
            && "Error: size mismatch!");
    x.resize(n);

#if SOLUTION
	MatrixXd A = buildA(a);
    x = A.fullPivLu().solve(b);
#else // TEMPLATE
    // TODO: solve system $A^{-1}*B$
#endif // TEMPLATE
}
/* SAM_LISTING_END_1 */

/* @brief Solve $Ax = b$.
 * Function with efficient implementation, using the solution of subproblem 4.
 * @param[in] a An $n$-dimensional vector to build the lower triangular matrix $A$
 * @param[in] b An $n$-dimensional vector
 * @param[out] x The $n$-dimensional vector $x = A^{-1}*b$
 */
/* SAM_LISTING_BEGIN_2 */
void solveA_fast(const VectorXd & a, const VectorXd & b, VectorXd & x) {
    // Size of b, which is the size of a
    int n = b.size();
    assert( n == a.size()
            && "Error: size mismatch!");
    x.resize(n);

#if SOLUTION
    x(0) = b(0)/a(0);
    for(int i=1; i<n; i++) {
		x(i) = (b(i) - b(i-1))/a(i);
	}
#else // TEMPLATE
    // TODO: solve system $A^{-1}*B$ efficiently
#endif // TEMPLATE
}
/* SAM_LISTING_END_2 */

/* @brief Test previous solutions
 */
int main() {
    unsigned int n = 9;
    // Compute with both solvers
    std::cout << "--> Check that the solvers are correct" << std::endl;
    VectorXd a = VectorXd::Random(n);
    VectorXd b = VectorXd::Random(n);
    VectorXd x1, x2, x;

    solveA(a,b,x1);
    solveA_fast(a,b,x2);

    std::cout << "Error = " << (x1 - x2).norm() << std::endl;

    // Compute runtimes of different solvers
    std::cout << "--> Runtime comparison of naive vs fast solver" << std::endl;
	// Number of repetitions
    unsigned int repeats = 3;


    // sizes will contain the size of the matrix
    // timings will contain the runtimes in seconds
    std::vector<double> sizes, timings, timings_eff;


    // Header
    std::cout << std::setw(20) << "n"
              << std::setw(20) << "time naive [s]"
              << std::setw(20) << "time fast [s]"
              << std::endl;

    // Loop over matrix size
    for(unsigned int k = 4; k <= 12; ++k) {
        // Timers
        Timer tm_naive, tm_fast;
        unsigned int n = pow(2,k);

        // Repeat test many times
        for(unsigned int r = 0; r < repeats; ++r) {
            a = VectorXd::Random(n);
            b = VectorXd::Random(n);

            // Compute runtime with naive solver
            tm_naive.start();
            solveA(a,b,x);
            tm_naive.stop();
            // Compute runtime with efficient solver
            tm_fast.start();
            solveA_fast(a,b,x);
            tm_fast.stop();
        }


        // Save results in a vector
        sizes.push_back(n); // save vector sizes
        timings.push_back(tm_naive.min());
        timings_eff.push_back(tm_fast.min());


        // Print runtimes
        std::cout << std::setw(20) << n
                  << std::scientific << std::setprecision(3)
                  << std::setw(20) << tm_naive.min()
                  << std::setw(20) << tm_fast.min()
                  << std::endl;
    }


    mgl::Figure fig;
    fig.title("Timings of naive solver");
    fig.ranges(2, 9000, 1e-8, 1e3);
    fig.setlog(true, true); // set loglog scale
    fig.plot(sizes, timings, " r+").label("original");
    fig.fplot("1e-9*x^3", "k|").label("O(n^3)");
    fig.xlabel("Vector size (n)");
    fig.ylabel("Time [s]");
    fig.legend(0, 1);
    fig.save("structuredLSE_timing.eps");
    fig.save("structuredLSE_timing.png");

    mgl::Figure fig2;
    fig2.title("Comparison of timings");
    fig2.ranges(2, 9000, 1e-8, 1e3);
    fig2.setlog(true, true); // set loglog scale
    fig2.plot(sizes, timings, " r+").label("original");
    fig2.plot(sizes, timings_eff, " b+").label("efficient");
    fig2.fplot("1e-9*x^3", "k|").label("O(n^3)");
    fig2.fplot("1e-7*x", "k-").label("O(n)");
    fig2.xlabel("Vector size (n)");
    fig2.ylabel("Time [s]");
    fig2.legend(0, 1);
    fig2.save("structuredLSE_comparison.eps");
    fig2.save("structuredLSE_comparison.png");

}
