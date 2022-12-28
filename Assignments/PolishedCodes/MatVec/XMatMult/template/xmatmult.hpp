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