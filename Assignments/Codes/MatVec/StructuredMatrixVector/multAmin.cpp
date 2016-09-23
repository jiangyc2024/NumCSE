#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"
#if INTERNAL
#include <figure/figure.hpp>
#endif // INTERNAL

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
#if SOLUTION
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
#else // TEMPLATE
    // TODO: Implement an efifcient version of multAmin
#endif // TEMPLATE
}
/* SAM_LISTING_END_3 */

int main(void) {
    // Testing correctness of the code
    unsigned int M = 10;
    VectorXd xa = VectorXd::Random(M);
    VectorXd ys, yf;

    multAmin(xa, yf);
    multAminSlow(xa, ys);
    // Error should be small
    std::cout << "||ys-yf|| = " << (ys - yf).norm() << std::endl;


#if SOLUTION
    /* SAM_LISTING_BEGIN_4 */
    // Timing from $2^4$ to $2^{13}$ repeating "nruns" times
    unsigned int nruns = 10;
#if INTERNAL
    std::vector<double> sizes, times_slow,
        times_slow_loops, times_fast;
#endif // INTERNAL

    std::cout << "--> Timings:" << std::endl;
    // Header, see iomanip documentation
    std::cout << std::setw(15)
              << "N"
              << std::scientific << std::setprecision(3)
              << std::setw(15) << "multAminSlown"
              << std::setw(15) << "multAminLoops"
              << std::setw(15) << "multAmin"
              << std::endl;
    // From $2^4$ to $2^{13}$
    for(unsigned int N = (1 << 4); N <= (1 << 13); N = N << 1) {
        Timer tm_slow, tm_slow_loops, tm_fast;
        // Compute runtime many times
        for(unsigned int r = 0; r < nruns; ++r) {
            VectorXd x = VectorXd::Random(N);
            VectorXd y;

            // Runtime of slow method
            tm_slow.start();
            multAminSlow(x, y);
            tm_slow.stop();

            // Runtime of slow method with loops
            tm_slow_loops.start();
            multAminLoops(x, y);
            tm_slow_loops.stop();

            // Runtime of fast method
            tm_fast.start();
            multAmin(x, y);
            tm_fast.stop();
        }

#if INTERNAL
        sizes.push_back(N);
        times_slow.push_back( tm_slow.min() );
        times_slow_loops.push_back( tm_slow_loops.min() );
        times_fast.push_back( tm_fast.min() );
#endif // INTERNAL

        std::cout << std::setw(15)
                  << N
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << tm_slow.min()
                  << std::setw(15) << tm_slow_loops.min()
                  << std::setw(15) << tm_fast.min()
                  << std::endl;
    }
    /* SAM_LISTING_END_4 */
#else // TEMPLATE
 // TODO: Time multAminSlow and multAmin
 // Repeat timings 10 times. Output times inseconds with
 // scientific notation and 3 digits
#endif //TEMPLATE

#if INTERNAL
  // Plotting
  mgl::Figure fig;
  fig.title("Timings of multAmin");
  fig.ranges(2, 9000, 1e-8, 1e3);
  fig.setlog(true, true); // set loglog scale
  fig.plot(sizes, times_slow, " r+")
      .label("multAmin (slow)");
  fig.plot(sizes, times_slow_loops, " b+")
      .label("multAmin (slow, with loops)");
  fig.plot(sizes, times_fast, " g+")
      .label("multAmin");
  fig.fplot("1e-9*x^2", "k|").label("O(n^2)");
  fig.fplot("1e-9*x", "k-").label("O(n)");
  fig.xlabel("Vector size (n)");
  fig.ylabel("Time [s]");
  fig.legend(0, 1);
  fig.save("multAmin_timing.eps");
  fig.save("multAmin_timing.png");
#endif // INTERNAL

#if SOLUTION
    // The following code is just for demonstration purposes.
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
    std::cout << "--> Test B = inv(A):" << std::endl;
    VectorXd x = VectorXd::Random(n), y;
    multAmin(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
    multAminSlow(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
    multAminLoops(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
#endif
}
