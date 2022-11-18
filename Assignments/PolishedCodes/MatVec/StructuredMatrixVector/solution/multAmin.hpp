#ifndef MULTAMIN_HPP
#define MULTAMIN_HPP
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "timer.h"

/**
 * @brief compute $\mathbf{A}\mathbf{x}$
 * $\mathbf{A}$ is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * @param x vector x for computation of A*x = y
 * @param y = A*x
 */
void multAminSlow(const Eigen::VectorXd &x, Eigen::VectorXd &y) {
  const unsigned int n = x.size();

  /* SAM_LISTING_BEGIN_1 */
  Eigen::VectorXd one = Eigen::VectorXd::Ones(n);
  Eigen::VectorXd linsp = Eigen::VectorXd::LinSpaced(n, 1, n);
  y = ((one * linsp.transpose()).cwiseMin(linsp * one.transpose())) * x;
  /* SAM_LISTING_END_1 */
}

/**
 * @brief compute $\mathbf{A}\mathbf{x}$
 * This function has optimal complexity.
 * $\mathbf{A}$ is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * @param x vector x for computation of A*x = y
 * @param y = A*x
 */
/* SAM_LISTING_BEGIN_2 */
void multAmin(const Eigen::VectorXd &x, Eigen::VectorXd &y) {
  const unsigned int n = x.size();
  y = Eigen::VectorXd::Zero(n);

  // TODO: (1-7.b) Fill in the entries of y.
  // Hint: Find an expression y(j) = w(j) + j*u(j) such that
  // the entries of w and u can be calculated recursively.
  // START
  Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(n);

  v(0) = x(n - 1);
  w(0) = x(0);

  for (unsigned int j = 1; j < n; ++j) {
    v(j) = v(j - 1) + x(n - j - 1);
    w(j) = w(j - 1) + (j + 1) * x(j);
  }
  for (unsigned int j = 0; j < n - 1; ++j) {
    y(j) = w(j) + v(n - j - 2) * (j + 1);
  }
  y(n - 1) = w(n - 1);
  // END
}
/* SAM_LISTING_END_2 */

void multAmin_runtime() {
  /* SAM_LISTING_BEGIN_3 */
  // Timing from $2^4$ to $2^{10}$ repeating "nruns" times
  constexpr unsigned int nruns = 10;

  std::cout << "--> Timings:" << std::endl;
  // Header, see iomanip documentation
  std::cout << std::setw(15) << "N" << std::scientific << std::setprecision(3)
            << std::setw(15) << "multAminSlow" << std::setw(15) << "multAmin"
            << std::endl;
  // From $2^4$ to $2^{10}$
  // Note: << and >> are the bitwise shift operators
  for (unsigned int N = (1 << 4); N <= (1 << 10); N = N << 1) {
    Timer tm_slow, tm_fast;
    // TODO: (1-7.d) Compute runtimes of multAminSlow(x,y) and
    // multAmin(x,y) with x = Eigen::VectorXd::Random(N). Repeat nruns times.
    // START
    for (unsigned int r = 0; r < nruns; ++r) {
      Eigen::VectorXd x = Eigen::VectorXd::Random(N);
      Eigen::VectorXd y;

      // Runtime of slow method
      tm_slow.start();
      multAminSlow(x, y);
      tm_slow.stop();

      // Runtime of fast method
      tm_fast.start();
      multAmin(x, y);
      tm_fast.stop();
    }
    // END

    std::cout << std::setw(15) << N << std::scientific << std::setprecision(3)
              << std::setw(15) << tm_slow.min() << std::setw(15)
              << tm_fast.min() << std::endl;
  }
  /* SAM_LISTING_END_3 */
}

/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXd multABunitv() {
  constexpr unsigned int n = 10;

  /* SAM_LISTING_BEGIN_5 */
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, n);
  for (unsigned int i = 0; i < n; ++i) {
    B(i, i) = 2;
    if (i < n - 1) B(i + 1, i) = -1.0;
    if (i > 0) B(i - 1, i) = -1.0;
  }
  B(n - 1, n - 1) = 1.0;
  /* SAM_LISTING_END_5 */

  Eigen::MatrixXd C(n, n);

  // TODO: (1-7.f) Set the columns of C to A*B*e_j for
  // j=0,...,n-1, and print C.
  // START
  Eigen::VectorXd y;
  // B*e_j is the j-th column of B*Id = B, so we need only
  // calculate A*B.col(j) for j=0,...,n-1.
  for (unsigned int i = 0; i < n; ++i) {
    multAmin(B.col(i), y);
    C.col(i) = y;
  }
  std::cout << "C = \n" << C << std::endl;
  // END

  return C;
}
/* SAM_LISTING_END_4 */

#endif
