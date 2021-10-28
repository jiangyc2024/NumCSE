////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <cmath>
#include <complex>
#include <iostream>

// \brief Compute complex root
// @param w complex number with non negative imaginary parts
// @param the square root of w with non negative real and imaginary parts

/* SAM_LISTING_BEGIN_0 */
std::complex<double> myroot(std::complex<double> w) {
  double x, y;
  double u = w.real();
  double v = w.imag();

  // TODO: (1-11.c) Compute the square root of w avoiding cancellation,
  // use only real arithmetic.

  // START

  // END

  return std::complex<double>(x, y);
}
/* SAM_LISTING_END_0 */
