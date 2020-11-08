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
MatrixXd toeplitz(const VectorXd &c, const VectorXd &r) {
  // Find size of Toeplitz matrix
  const int m = c.size();
  const int n = r.size();
  MatrixXd T(m, n);
  // TODO: build Toeplitz matrix $\VT$
  // START
  for (int i = 0; i < n; ++i) {
    T.col(i).tail(m - i) = c.head(m - i);
  }
  for (int i = 0; i < m; ++i) {
    T.row(i).tail(n - i - 1) = r.segment(1, n - i - 1);
  }  // Do not reassign the diagonal!
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
VectorXd toepmatmult(const VectorXd &c, const VectorXd &r, const VectorXd &x) {
  assert(c.size() == r.size() && c.size() == x.size() &&
         "c, r, x have different lengths!");
  return toeplitz(c, r) * x;
}
/* SAM_LISTING_END_0 */

/* @brief Do something...
 * @param[in] c An $n$-dimensional vector
 * @param[in] r An $n$-dimensional vector
 * @param[in] x An $n$-dimensional vector
 * @param[out] y An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd toepmult(const VectorXd &c, const VectorXd &r, const VectorXd &x) {
  assert(c.size() == r.size() && c.size() == x.size() &&
         "c, r, x have different lengths!");
  int n = c.size();
  // Complex arithmetic required here
  VectorXcd cr_tmp = c.cast<std::complex<double>>();
  VectorXcd x_tmp = x.cast<std::complex<double>>();
  // Prepare vector encoding circulant matrix $\VC$
  cr_tmp.conservativeResize(2 * n);
  cr_tmp.tail(n) = VectorXcd::Zero(n);
  cr_tmp.tail(n - 1).real() = r.tail(n - 1).reverse();
  // Zero padding
  x_tmp.conservativeResize(2 * n);
  x_tmp.tail(n) = VectorXcd::Zero(n);
  // Periodic discrete convolution from \lref{cpp:pconvfft}
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
VectorXd ttmatsolve(const VectorXd &h, const VectorXd &y) {
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
VectorXd ttrecsolve(const VectorXd &h, const VectorXd &y, int l) {
  assert(h.size() == y.size() && "h and y have different lengths!");
  // Result vector
  VectorXd x;
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
VectorXd ttsolve(const VectorXd &h, const VectorXd &y) {
  assert(h.size() == y.size() && "h and y have different lengths!");
  int n = h.size();

  VectorXd x;

  // TODO: wrap 'ttrecsolve' for any size $n$
  // START
  int l = std::ceil(std::log(n));
  int m = std::pow(2, l);

  VectorXd h_tmp = h;
  h_tmp.conservativeResize(m);

  VectorXd y_tmp = y;
  y_tmp.conservativeResize(m);

  x = ttrecsolve(h_tmp, y_tmp, l);
  x.conservativeResize(n);
  x.bottomRows(n - m).setZero();
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

  for (unsigned int l = 3; l <= 9; l += 1) {
    // vector size
    unsigned int n = std::pow(2, l);
    // save vector size n
    vec_size.push_back(n);

    // number of repetitions
    unsigned int repeats = 3;
    // TODO: (5-4.g) Perform a runtime comparison by repeating time computation
    // 'repeats' times
    // START
    Timer tm_matmult, tm_mult, tm_ttmat, tm_ttrec;
    // repeat test 'repeats' times
    for (unsigned int rr = 0; rr < repeats; ++rr) {
      // create runtime test using randome vectors and given vector h
      VectorXd h = VectorXd::LinSpaced(n, 1, n).cwiseInverse();
      VectorXd c = VectorXd::Random(n);
      VectorXd r = VectorXd::Random(n);
      VectorXd x = VectorXd::Random(n);
      VectorXd y = VectorXd::Random(n);
      r(0) = c(0);

      // compute times for toepmatmult implementation
      tm_matmult.start();
      toepmatmult(c, r, x);
      tm_matmult.stop();
      // compute times for toepmult implementation
      tm_mult.start();
      toepmult(c, r, x);
      tm_mult.stop();

      // compute times for ttmatsolve implementation
      tm_ttmat.start();
      ttmatsolve(h, y);
      tm_ttmat.stop();
      // compute times
      tm_ttrec.start();
      ttrecsolve(h, y, l);
      tm_ttrec.stop();
    }
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

    /* DO NOT CHANGE */
    // create plot
    plot(vec_size, elap_time_mult, elap_time_matmult, "./cx_out/fig1.png",
         "toepmult", "toepmatmult");
    plot(vec_size, elap_time_ttrec, elap_time_ttmat, "./cx_out/fig2.png",
         "ttrecsolve", "ttmatsolve");
  }
}
/* SAM_LISTING_END_6 */

// Additional

// rename long variable name to duration_t (easy to change)
using duration_t = std::chrono::nanoseconds;

/* \brief Compute runtime of $F$, repeating test "repeats" times
 * Will return minimal runtime.
 * This function uses "crhono".
 * \tparam Function type of F, must have an operator()
 * \param[in] F Function for which you want to measure runtime.
 * \param[in] repeats Number of repetitions.
 */
template <class Function>
duration_t timing(const Function &F, int repeats = 10) {
  // Shortcut for time_point
  using time_point_t = std::chrono::high_resolution_clock::time_point;

  // Loop many times
  duration_t min_elapsed;
  for (int r = 0; r < repeats; r++) {
    // Start clock (MATLAB: tic)
    time_point_t start = std::chrono::high_resolution_clock::now();

    // Run function
    F();

    // Stop clock (MATLAB: toc) and measure difference
    duration_t elapsed = std::chrono::duration_cast<duration_t>(
        std::chrono::high_resolution_clock::now() - start);

    // Compute min between all runs
    min_elapsed = r == 0 ? elapsed : std::min(elapsed, min_elapsed);
  }

  return min_elapsed;
}

/* \brief Compute timing using chrono
 * Also demonstrate use of lambda functions
 */
void runtime_toeplitz_with_chrono() {
  // table header
  std::cout << std::setw(8) << "n" << std::setw(15) << "toepmatmult"
            << std::setw(15) << "toepmult" << std::setw(20) << "ttmatsolve"
            << std::setw(15) << "ttrecsolve" << std::endl;

  // vector size
  unsigned int n;
  // repeat test 'repeats' times
  for (unsigned int l = 3; l <= 9; ++l) {
    // vectore size
    n = pow(2, l);

    // create runtime test using randome vectors and given vector h
    VectorXd h = VectorXd::LinSpaced(n, 1, n).cwiseInverse();
    VectorXd c = VectorXd::Random(n);
    VectorXd r = VectorXd::Random(n);
    VectorXd x = VectorXd::Random(n);
    VectorXd y = VectorXd::Random(n);
    r(0) = c(0);

    // Call "timing", using a lambda function for F
    // Remember: we cannot pass arrow\_matrix\_2\_times\_x directly to timing
    // the timing function expects a n object with operator()(void)
    duration_t elap_time_matmult =
        timing([&c, &r, &x]() { toepmatmult(c, r, x); }, 3);
    duration_t elap_time_mult =
        timing([&c, &r, &x]() { toepmult(c, r, x); }, 3);
    duration_t elap_time_ttmat = timing([&y, &h]() { ttmatsolve(y, h); }, 3);
    duration_t elap_time_ttrec =
        timing([&y, &h, &l]() { ttrecsolve(y, h, l); }, 3);

    // output timings chrono excercice
    // print the results: toepmult vs toepmatmult
    std::cout << std::setw(8) << n << std::scientific << std::setprecision(3)
              << std::setw(15) << elap_time_matmult.count() * 1e-9  // ns to s
              << std::setw(15) << elap_time_mult.count() * 1e-9     // ns to s
              << std::setw(15) << elap_time_ttmat.count() * 1e-9    // ns to s
              << std::setw(15) << elap_time_ttrec.count() * 1e-9    // ns to s
              << std::endl;
  }
}

#endif
