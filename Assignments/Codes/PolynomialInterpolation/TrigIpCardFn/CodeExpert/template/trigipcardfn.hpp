#ifndef TRIGIPCARDFN_HPP
#define TRIGIPCARDFN_HPP

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iostream>

#include <unsupported/Eigen/FFT>

using namespace Eigen;

/*
 * @brief Efficient FFT-based computation of coefficients in expansion
 *  \eqref{eq:trigpreal} for a trigonometric interpolation polynomial in
 * equidistant points \Blue{$(\frac{j}{2n+1},y_j)$}, \Blue{$j=0,\ldots,2n$}. IN
 *  : \texttt{y} has to be a row vector of odd length, return values are column
 *  vectors
 *      \texttt{a}, \texttt{b} will be used to save the expansion coefficients
 */

void trigipequid(const VectorXd& y, VectorXcd& a, VectorXcd& b) {
  const unsigned N = y.size();
  if (N % 2 != 1) {
    std::cerr << "Number of points must be odd!\n";
    return;
  }
  const unsigned n = (N - 1) / 2;
  // prepare data for fft
  std::complex<double> i(0, 1);  // imaginary unit
  std::vector<std::complex<double> > f(N);
  std::vector<std::complex<double> > c;

  for (unsigned k = 0; k < N; ++k) {
    // see \eqref{tip:FM}
    f[k] = y(k) * std::exp(2 * M_PI * i * (double(n) / N * k));
  }
  FFT<double> fft;
  fft.fwd(c, f);  // -> c = fft(f);

  // From \eqref{eq:trigpcomp}: \Blue{$\alpha_j =
  // \frac{1}{2}(\gamma_{n-j}+\gamma_{n+j})$} and \Blue{$\beta_j =
  // \frac{1}{2i}(\gamma_{n-j}-\gamma_{n+j})$}, \Blue{$j=1,\ldots,n$},
  // \Blue{$\alpha_0 = \gamma_n$}
  a = VectorXcd(n + 1);
  b = VectorXcd(n);

  a(0) = c[n];
  for (unsigned l = 1; l <= n; ++l) {
    a(l) = c[n - l] + c[n + l];
    b(l - 1) = -i * (c[n - l] - c[n + l]);
  }
  // dont forget scaling factor of forward FFT!
  a /= N;
  b /= N;
}

/*
* @brief Evaluation of trigonometric interpolation polynomial through
* \Blue{$(\frac{j}{2n+1},y_j)$}, \Blue{$j=0,\ldots,2n$} in equidistant points
* \Blue{$\frac{k}{N}$}, \Blue{$k=0,N-1$} IN : \texttt{y} = vector of values to
* be interpolated
*      \texttt{q} (COMPLEX!) will be used to save the return values
*/
void trigpolyvalequid(const VectorXd y, const int M, VectorXd& q) {
  const int N = y.size();
  if (N % 2 == 0) {
    std::cerr << "Number of points must be odd!\n";
    return;
  }
  const int n = (N - 1) / 2;
  // computing coefficient \Blue{$\gamma_j$}, see \eqref{tip:FM}
  VectorXcd a, b;
  trigipequid(y, a, b);

  std::complex<double> i(0, 1);
  VectorXcd gamma(2 * n + 1);
  gamma(n) = a(0);
  for (int k = 0; k < n; ++k) {
    gamma(k) = 0.5 * (a(n - k) + i * b(n - k - 1));
    gamma(n + k + 1) = 0.5 * (a(k + 1) - i * b(k));
  }

  // zero padding
  VectorXcd ch(M);
  ch << gamma, VectorXcd::Zero(M - (2 * n + 1));

  // build conjugate fourier matrix
  FFT<double> fft;
  const VectorXcd chCon = ch.conjugate();
  const VectorXcd v = fft.fwd(chCon).conjugate();

  // multiplicate with conjugate fourier matrix
  VectorXcd q_complex = VectorXcd(M);
  for (int k = 0; k < M; ++k) {
    q_complex(k) = v(k) * std::exp(-2. * k * n * M_PI / M * i);
  }
  // complex part is zero up to machine precision, cut off!
  q = q_complex.real();
}

/*!
 * @brief trigIpL Compute $\lambda(n)$.
 *
 * @param[in] n $2*n+1$ will be the number of basis polynomials.
 * @param[out] Value $\lambda(n)$.
 */
/* SAM_LISTING_BEGIN_1 */
double trigIpL(std::size_t n) {
  
  // TO DO: write a function that approximatly computes the Lebesgue constant $\lambda(n)$ for n = 2$^k$, k = 2,3, ..6
  // START

  // END
}
/* SAM_LISTING_END_1 */

#endif
