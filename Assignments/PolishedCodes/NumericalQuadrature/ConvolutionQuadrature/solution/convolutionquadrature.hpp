#ifndef CONVOLUTIONQUADRATURE_HPP
#define CONVOLUTIONQUADRATURE_HPP
/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <unsupported/Eigen/FFT>
#include <vector>

// DFT-based computation of convolution quadrature weights.
/* SAM_LISTING_BEGIN_1 */
template <typename FFUNCTOR>
Eigen::VectorXcd compute_cq_weights(FFUNCTOR&& F, unsigned int N, double tau,
                                    double r = 0.0) {
  using Comp = std::complex<double>;
  assert(N >= 1 && "At least two quadrature points have to be used!");
  // Default choice for circle radius for integration contour
  if (r == 0.0) {
    r = std::pow(std::numeric_limits<double>::epsilon(), 0.5 / (N + 1));
  }
  Eigen::VectorXcd w(N + 1);
  const Comp iu = 2.0 * Comp{0.0, 1.0} * M_PI;  // 2*pi*i
  // TODO: (7-15.c) Compute the approximate convolution quadrature weights by
  // means of the equidistant trapezoidal rule.
  // START
  Eigen::VectorXcd v(N + 1);
  // Sampling of F on a circle in the complex plane, see \prbeqref{eq:x}
  const Comp expfac = iu / (static_cast<double>(N + 1));
  for (unsigned int k = 0; k <= N; ++k) {
    const Comp z =
        (Comp(1.0, 0.0) - r * std::exp(expfac * static_cast<double>(k))) / tau;
    v[k] = F(z);
  }
  // Compute DFT of the vector v
  Eigen::FFT<double> fft;
  w = fft.fwd(v);  // DFT in \eigen
  // Finally perform the scaling of the components
  double fac = N + 1;
  for (unsigned int k = 0; k <= N; ++k) {
    w[k] /= fac;
    fac *= r;
  }
  // END
  return w;
}
/* SAM_LISTING_END_1 */

// Exact convolution quadrature weights for F the square root function
Eigen::VectorXd compute_cq_weights_sqrt(unsigned int N, double tau) {
  Eigen::VectorXd w(N + 1);
  double fac = 1.0 / std::sqrt(tau);
  for (unsigned int j = 0; j <= N; ++j) {
    w[j] = fac;
    fac *= (-(0.5 - j) / (j + 1.0));
  }
  return w;
}

#endif