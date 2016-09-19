#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* \brief compute $\mathbf{A}\mathbf{x}$
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminSlow(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    /* SAM_LISTING_BEGIN_1 */
    VectorXd one = VectorXd::Ones(n);
    VectorXd linsp = VectorXd::LinSpaced(n,1,n);
    y = ( ( one * linsp.transpose() )
          .cwiseMin( linsp * one.transpose()) ) * x;
    /* SAM_LISTING_END_1 */
}

/* \brief compute $\mathbf{A}\mathbf{x}$
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * Instead of a "Matlab style" construcion of the product,
 * we use simple loops.
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminLoops(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    MatrixXd A(n,n);

    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            A(i,j) = std::min(i+1,j+1);
        }
    }
    y = A * x;
}

/* \brief compute $\mathbf{A}\mathbf{x}$
 * This function has optimal complexity.
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
/* SAM_LISTING_BEGIN_3 */
void multAmin(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();
    y = VectorXd::Zero(n);
    VectorXd v = VectorXd::Zero(n);
    VectorXd w = VectorXd::Zero(n);

    v(0) = x(n-1);
    w(0) = x(0);

    for(unsigned int j = 1; j < n; ++j) {
        v(j) = v(j-1) + x(n-j-1);
        w(j) = w(j-1) + (j+1)*x(j);
    }
    for(unsigned int j = 0; j < n-1; ++j) {
        y(j) = w(j) + v(n-j-2)*(j+1);
    }
    y(n-1) = w(n-1);
}
/* SAM_LISTING_END_3 */

int main(void) {
    // Timing from $2^4$ to $2^13$ repeating "nruns" times
    unsigned int nruns = 10;
    std::vector<double> sizes, times_slow,
        times_slow_loops, times_fast;
    for(unsigned int N = (1 << 4); N <= (1 << 13); N = N << 1) {
        Timer tm_slow, tm_slow_loops, tm_fast;
        for(unsigned int r = 0; r < nruns; ++r) {
            VectorXd x = VectorXd::Random(N);
            VectorXd y;

            tm_slow.start();
            multAminSlow(x, y);
            tm_slow.stop();

            tm_slow_loops.start();
            multAminLoops(x, y);
            tm_slow_loops.stop();

            tm_fast.start();
            multAmin(x, y);
            tm_fast.stop();
        }

        sizes.push_back(N);
        times_slow.push_back( tm_slow.min() );
        times_slow_loops.push_back( tm_slow_loops.min() );
        times_fast.push_back( tm_fast.min() );

        std::cout << std::setw(15)
                  << N
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << tm_slow.min()
                  << std::setw(15) << tm_slow_loops.min()
                  << std::setw(15) << tm_fast.min()
                  << std::endl;
    }


    // The following code is kust for demonstration purposes.
    // Build Matrix B with dimension 10x10 such that B = inv(A)
    unsigned int n = 10;
    /* SAM_LISTING_BEGIN_2 */
    MatrixXd B = MatrixXd::Zero(n,n);
    for(unsigned int i = 0; i < n; ++i) {
        B(i,i) = 2;
        if(i < n-1) B(i+1,i) = -1;
        if(i > 0) B(i-1,i) = -1;
    }
    B(n-1,n-1) = 1;
    /* SAM_LISTING_END_2 */
    std::cout << "B = " << std::endl
              << B << std::endl;

    // Check that B = inv(A) (up to machine precision)
    VectorXd x = VectorXd::Random(n), y;
    multAmin(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
    multAminSlow(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
    multAminLoops(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
}
