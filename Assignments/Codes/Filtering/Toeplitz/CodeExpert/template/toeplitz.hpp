#ifndef TOEPLITZ_HPP
#define TOEPLITZ_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>

#include "matplotlibcpp.h"
#include "pconvfft.hpp"
#include "plot.hpp"
#include "timer.h"

using namespace Eigen;

/* @brief Build a Toeplitz matrix $\VT$ from $\Vc$ and $\Vr$
 * @param[in] c An $m$-dimensional vector, first column of $\VT$
 * @param[in] r An $n$-dimensional vector, first row of $\VT$
 * @param[out] T The $m \times n$ Toeplitz matrix from $\Vc$ and $\Vr$
 */
/* SAM_LISTING_BEGIN_5 */
MatrixXd toeplitz(const VectorXd& c, const VectorXd& r) {
  if (c(0) != r(0)) {
    std::cerr << "First entries of c and r are different!" << std::endl
              << "We assign the first entry of c to the diagonal" << std::endl;
  }

  // Initialization
  int m = c.size();
  int n = r.size();
  MatrixXd T(m, n);

  // TODO: build Toeplitz matrix $\VT$
  // START

  // END

  return T;
}
/* SAM_LISTING_END_5 */

/* @brief Do something...
 * @param[in] c An $n$-dimensional vector
 * @param[in] r An $n$-dimensional vector
 * @param[in] x An $n$-dimensional vector
 * @param[out] y An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd toepmatmult(const VectorXd& c, const VectorXd& r, const VectorXd& x) {
  assert(c.size() == r.size() && c.size() == x.size() &&
         "c, r, x have different lengths!");

  MatrixXd T = toeplitz(c, r);

  VectorXd y = T * x;

  return y;
}
/* SAM_LISTING_END_0 */

/* @brief Do something...
 * @param[in] c An $n$-dimensional vector
 * @param[in] r An $n$-dimensional vector
 * @param[in] x An $n$-dimensional vector
 * @param[out] y An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd toepmult(const VectorXd& c, const VectorXd& r, const VectorXd& x) {
  assert(c.size() == r.size() && c.size() == x.size() &&
         "c, r, x have different lengths!");
  int n = c.size();

  VectorXcd cr_tmp = c.cast<std::complex<double>>();
  cr_tmp.conservativeResize(2 * n);
  cr_tmp.tail(n) = VectorXcd::Zero(n);
  cr_tmp.tail(n - 1).real() = r.tail(n - 1).reverse();

  VectorXcd x_tmp = x.cast<std::complex<double>>();
  x_tmp.conservativeResize(2 * n);
  x_tmp.tail(n) = VectorXcd::Zero(n);

  VectorXd y = pconvfft(cr_tmp, x_tmp).real();
  y.conservativeResize(n);

  return y;
}
/* SAM_LISTING_END_1 */

/* @brief Do something...
 * @param[in] h An $n$-dimensional vector
 * @param[in] y An $n$-dimensional vector
 * @param[out] x An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd ttmatsolve(const VectorXd& h, const VectorXd& y) {
  assert(h.size() == y.size() && "h and y have different lengths!");
  int n = h.size();

  VectorXd h_tmp = VectorXd::Zero(n);
  h_tmp(0) = h(0);

  MatrixXd T = toeplitz(h, h_tmp);

  VectorXd x = T.fullPivLu().solve(y);

  return x;
}
/* SAM_LISTING_END_2 */

/* @brief Do something...
 * @param[in] h An $n$-dimensional vector
 * @param[in] y An $n$-dimensional vector
 * @param[in] l An integer
 * @param[out] x An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_3 */
VectorXd ttrecsolve(const VectorXd& h, const VectorXd& y, int l) {
  assert(h.size() == y.size() && "h and y have different lengths!");

  VectorXd x;

  if (l == 0) {
    x.resize(1);
    x(0) = y(0) / h(0);
  } else {
    int n = std::pow(2, l);
    int m = n / 2;

    assert(h.size() == n && y.size() == n &&
           "h and y have length different from 2^l!");

    VectorXd x1 = ttrecsolve(h.head(m), y.head(m), l - 1);
    VectorXd y2 = y.segment(m, m) -
                  toepmult(h.segment(m, m), h.segment(1, m).reverse(), x1);
    VectorXd x2 = ttrecsolve(h.head(m), y2, l - 1);

    x.resize(n);
    x.head(m) = x1;
    x.tail(m) = x2;
  }

  return x;
}
/* SAM_LISTING_END_3 */

/* @brief Wrapper for 'ttrecsolve' for any size $n$
 * @param[in] h An $n$-dimensional vector
 * @param[in] y An $n$-dimensional vector
 * @param[out] x An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_4 */
VectorXd ttsolve(const VectorXd& h, const VectorXd& y) {
  assert(h.size() == y.size() && "h and y have different lengths!");
  int n = h.size();

  VectorXd x;

  // TODO: wrap 'ttrecsolve' for any size $n$
  // START

  // END

  return x;
}
/* SAM_LISTING_END_4 */

/* \brief Compute the runtime comparison of
 * toepmatmult vs toepmult and ttmatsolve vs ttrecsolve
 * Repeat tests 10 times, and output the minimal runtime amongst all times.
 */
/* SAM_LISTING_BEGIN_6 */
void runtime_toeplitz() {
  // memory allocation for plot
  std::vector<double> vec_size;
  std::vector<double> elap_time_matmult, elap_time_mult, elap_time_ttmat,
      elap_time_ttrec;

  // header for the results to print out
  std::cout << std::setw(8) << "n" << std::setw(15) << "toepmatmult"
            << std::setw(15) << "toepmult" << std::setw(20) << "ttmatsolve"
            << std::setw(15) << "ttrecsolve" << std::endl;

  // vector size
  unsigned int n;
  Timer tm_matmult, tm_mult, tm_ttmat, tm_ttrec;
  for (unsigned int l = 3; l <= 8; l += 1) {
    // vector size
    n = std::pow(2, l);
    // save vector size n
    vec_size.push_back(n);

    // number of repetitions
    unsigned int repeats = 3;

    // TODO: (5-4.g)  Perform a runtme comparison by repeating time computation
    // 'repeats' times 
    // START

    // END
  }

  // print the results: toepmult vs toepmatmult
  std::cout << std::setw(8) << n << std::scientific << std::setprecision(3)
            << std::setw(15) << tm_mult.min() << std::setw(15)
            << tm_matmult.min() << std::setw(20) << tm_ttmat.min()
            << std::setw(15) << tm_ttrec.min() << std::endl;

  // save elapsed time for plot: toepmatmult vs toepmult
  elap_time_matmult.push_back(tm_matmult.min());
  elap_time_mult.push_back(tm_mult.min());
  // save elapsed time for plot: ttmatsove vs ttrecsolve
  elap_time_ttmat.push_back(tm_matmult.min());
  elap_time_ttrec.push_back(tm_mult.min());
  /* DO NOT CHANGE */
  // create plot
  plot(vec_size, elap_time_mult, elap_time_matmult, "./cx_out/fig1.png",
       "toepmult", "toepmatmult");
  plot(vec_size, elap_time_ttrec, elap_time_ttmat, "./cx_out/fig2.png",
       "ttrecsolve", "ttmatsolve");
}
/* SAM_LISTING_END_6 */

#endif
