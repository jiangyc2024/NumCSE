#ifndef XMATMULT_HPP
#define XMATMULT_HPP

#include <Eigen/Dense>
#include <cassert>
#include <iomanip>
#include <iostream>

#include "timer.h"

/**
 * \brief Efficiently compute the Matrix-Vector-Product
 * $x=Ay$, where A(i,j) is zero except for values on the
 * diagonal and antidiagonal. There the values are given
 * by the entries of the vector a.
 *
 * \param a vector storing the elements of $A$
 * \param y vector to multiply $A$ with
 * \param x result vector of $Ay=x$
 */
/* SAM_LISTING_BEGIN_1 */
template <class Vector>
void xmatmult(const Vector& a, const Vector& y, Vector& x) {
  assert(a.size() == y.size() && a.size() == x.size() &&
         "Input vector dimensions must match");
  const unsigned int n = a.size();

  // TODO: (1-13.a) Efficiently compute the matrix-vector product for a matrix
  // with nonzeros on diagonal and antidiagonal.
  // START

  // looping over half of a,
  // if n is odd we treat the middle element differently
  for (unsigned int i = 0; i < n / 2; ++i) {  // integer division
    x(i) = a(i) * y(i) + a(n - i - 1) * y(n - i - 1);
    x(n - i - 1) = x(i);
  }

  if (n % 2) {  // n is odd
    x(n / 2) = a(n / 2) * y(n / 2);
  }
  // END
}
/* SAM_LISTING_END_1 */

/**
 * \brief Compares the runtimes of the efficient
 * multipication to the normal Matrix-Vector product
 */
/* SAM_LISTING_BEGIN_2 */
void compare_times() {
  std::cout << "Measuring runtimes for comparison" << std::endl;
  constexpr unsigned int repeats = 3;
  Timer t_fast, t_slow;
  Eigen::MatrixXd results(10, 3);
  // TODO: (1-13.c) Tabulate the runtimes of the efficient mvmult in comparison
  // to the normal one.
  // START
  for (unsigned int k = 14; k > 4; k--) {
    // bitshift operator '<<': 1<<3 == pow(2,3)
    const unsigned int n = 1 << k;
    Eigen::VectorXd a, y, x(n);
    a = y = Eigen::VectorXd::Random(n);

    Eigen::MatrixXd A = a.asDiagonal();
    for (unsigned int i = 0; i < n; ++i) {
      A(n - i - 1, i) = A(i, i);
    }

    // measure multiple times
    for (unsigned int i = 0; i < repeats; i++) {
      t_fast.start();
      xmatmult(a, y, x);
      t_fast.stop();

      t_slow.start();
      x = A * y;
      t_slow.stop();
    }
    results(k - 5, 0) = n;
    results(k - 5, 1) = t_slow.min();
    results(k - 5, 2) = t_fast.min();
  }
  // END

  // print results
  std::cout << std::setw(8) << "n" << std::setw(15) << "original"
            << std::setw(15) << "efficient" << std::setprecision(5)
            << std::endl;
  for (unsigned int i = 0; i < results.rows(); i++) {
    std::cout << std::setw(8) << results(i, 0) << std::setw(15) << results(i, 1)
              << " s" << std::setw(15) << results(i, 2) << " s" << std::endl;
  }
}
/* SAM_LISTING_END_2 */

#endif