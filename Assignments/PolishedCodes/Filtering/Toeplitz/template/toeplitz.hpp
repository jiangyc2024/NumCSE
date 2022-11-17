#ifndef TOEPLITZ_HPP
#define TOEPLITZ_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>

#include "pconvfft.hpp"
#include "plot.hpp"
#include "timer.h"

/**
 * @brief Build a Toeplitz matrix $\VT$ from $\Vc$ and $\Vr$
 *
 * @param c An $m$-dimensional vector, first column of $\VT$
 * @param r An $n$-dimensional vector, first row of $\VT$
 * @return Eigen::MatrixXd The $m \times n$ Toeplitz matrix from $\Vc$ and $\Vr$
 */
/* SAM_LISTING_BEGIN_5 */
Eigen::MatrixXd toeplitz(const Eigen::VectorXd &c, const Eigen::VectorXd &r) {
  // Find size of Toeplitz matrix
  const unsigned int m = c.size();
  const unsigned int n = r.size();
  Eigen::MatrixXd T = Eigen::MatrixXd::Zero(m, n);

  // TODO: (4-4.b) Build Toeplitz matrix $\VT$
  // START

  // END

  return T;
}
/* SAM_LISTING_END_5 */

/**
 * @brief Do something...
 *
 * @param c An $n$-dimensional vector
 * @param r An $n$-dimensional vector
 * @param x An $n$-dimensional vector
 * @return Eigen::VectorXd An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd toepmatmult(const Eigen::VectorXd &c, const Eigen::VectorXd &r,
                            const Eigen::VectorXd &x) {
  assert(c.size() == r.size() && c.size() == x.size() &&
         "c, r, x have different lengths!");
  return toeplitz(c, r) * x;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Do something...
 *
 * @param c An $n$-dimensional vector
 * @param r An $n$-dimensional vector
 * @param x An $n$-dimensional vector
 * @return Eigen::VectorXd An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd toepmult(const Eigen::VectorXd &c, const Eigen::VectorXd &r,
                         const Eigen::VectorXd &x) {
  assert(c.size() == r.size() && c.size() == x.size() &&
         "c, r, x have different lengths!");
  const unsigned int n = c.size();

  // Complex arithmetic required here
  Eigen::VectorXcd cr_tmp = c.cast<std::complex<double>>();
  Eigen::VectorXcd x_tmp = x.cast<std::complex<double>>();

  // Prepare vector encoding circulant matrix $\VC$
  cr_tmp.conservativeResize(2 * n);
  cr_tmp.tail(n) = Eigen::VectorXcd::Zero(n);
  cr_tmp.tail(n - 1).real() = r.tail(n - 1).reverse();

  // Zero padding
  x_tmp.conservativeResize(2 * n);
  x_tmp.tail(n) = Eigen::VectorXcd::Zero(n);

  // Periodic discrete convolution from \lref{cpp:pconvfft}
  Eigen::VectorXd y = pconvfft(cr_tmp, x_tmp).real();
  y.conservativeResize(n);

  return y;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Do something...
 *
 * @param h An $n$-dimensional vector
 * @param y An $n$-dimensional vector
 * @return Eigen::VectorXd An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd ttmatsolve(const Eigen::VectorXd &h, const Eigen::VectorXd &y) {
  assert(h.size() == y.size() && "h and y have different lengths!");
  const unsigned int n = h.size();

  Eigen::VectorXd h_tmp = Eigen::VectorXd::Zero(n);
  h_tmp(0) = h(0);

  Eigen::MatrixXd T = toeplitz(h, h_tmp);

  Eigen::VectorXd x = T.fullPivLu().solve(y);

  return x;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Do something...
 *
 * @param h An $n$-dimensional vector
 * @param y An $n$-dimensional vector
 * @param l An integer
 * @return Eigen::VectorXd An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd ttrecsolve(const Eigen::VectorXd &h, const Eigen::VectorXd &y,
                           int l) {
  assert(h.size() == y.size() && "h and y have different lengths!");
  // Result vector
  Eigen::VectorXd x;
  // Trivial case of asn 1x1 LSE
  if (l == 0) {
    x.resize(1);
    x(0) = y(0) / h(0);
  } else {
    int n = std::pow(2, l);
    int m = n / 2;

    // Check matching length of vectors
    assert(h.size() == n && y.size() == n &&
           "h and y have length different from 2^l!");

    Eigen::VectorXd x1 = ttrecsolve(h.head(m), y.head(m), l - 1);
    Eigen::VectorXd y2 =
        y.segment(m, m) -
        toepmult(h.segment(m, m), h.segment(1, m).reverse(), x1);
    Eigen::VectorXd x2 = ttrecsolve(h.head(m), y2, l - 1);

    x.resize(n);
    x.head(m) = x1;
    x.tail(m) = x2;
  }

  return x;
}
/* SAM_LISTING_END_3 */

/**
 * @brief Wrapper for 'ttrecsolve' for any size $n$
 *
 * @param h An $n$-dimensional vector
 * @param y An $n$-dimensional vector
 * @return Eigen::VectorXd An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd ttsolve(const Eigen::VectorXd &h, const Eigen::VectorXd &y) {
  assert(h.size() == y.size() && "h and y have different lengths!");
  const unsigned int n = h.size();

  Eigen::VectorXd x;

  // TODO: (4-4.i) Wrap 'ttrecsolve' for any size $n$
  // START

  // END

  return x;
}
/* SAM_LISTING_END_4 */

/**
 * @brief Compute the runtime comparison of
 * toepmatmult vs toepmult and ttmatsolve vs ttrecsolve
 * Repeat tests 10 times, and output the minimal runtime amongst all times.
 *
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

  for (unsigned int l = 3; l <= 9; l += 1) {
    // vector size
    unsigned int n = std::pow(2, l);
    // save vector size n
    vec_size.push_back(n);

    // number of repetitions
    constexpr unsigned int repeats = 3;
    Timer tm_matmult, tm_mult, tm_ttmat, tm_ttrec;
    // TODO: (4-4.g) Perform a runtime comparison by repeating time computation
    // 'repeats' times
    // START

    // END

    // print the results: toepmult vs toepmatmult
    std::cout << std::setw(8) << n << std::scientific << std::setprecision(3)
              << std::setw(15) << tm_matmult.min() << std::setw(15)
              << tm_mult.min() << std::setw(20) << tm_ttmat.min()
              << std::setw(15) << tm_ttrec.min() << std::endl;

    // save elapsed time for plot: toepmatmult vs toepmult
    elap_time_matmult.push_back(tm_matmult.min());
    elap_time_mult.push_back(tm_mult.min());
    // save elapsed time for plot: ttmatsove vs ttrecsolve
    elap_time_ttmat.push_back(tm_ttmat.min());
    elap_time_ttrec.push_back(tm_ttrec.min());
  }

  // create plot
  plot(vec_size, elap_time_mult, elap_time_matmult, "./cx_out/fig1.png",
       "toepmult", "toepmatmult");
  plot(vec_size, elap_time_ttrec, elap_time_ttmat, "./cx_out/fig2.png",
       "ttrecsolve", "ttmatsolve");
}
/* SAM_LISTING_END_6 */

#endif
