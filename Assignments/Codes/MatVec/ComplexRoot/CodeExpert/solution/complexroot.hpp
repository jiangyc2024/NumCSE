//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <complex>
#include <cmath>

//! \brief Compute complex root
//! \param[in] w complex number with non negative imaginary parts
//! \return the square root of w with non negative real and imaginary parts

/* SAM_LISTING_BEGIN_0 */
std::complex<double> myroot( std::complex<double> w ) {
    double x,y;
    double u = w.real();
    double v = w.imag();
    
    // TODO: (2-11.c) Compute the square root of w avoiding cancellation,
    // use only real arithmetic.
    // START
    // w is actually a real number, we are done
    if (v==0) return std::sqrt(u);

    if (u > 0) {
        // Compute the real part, since stable
        x = std::sqrt((std::sqrt(u*u+v*v)+u)/2.);
        // Use product to compute imaginary part
        y = v/(2*x);
    } else {
        // Compute the imaginary part, since stable
        y = std::sqrt((std::sqrt(u*u+v*v)-u)/2.);
        // Use product to compute real part
        x = v/(2*y);
    }
    // END 
    
    return std::complex<double> (x,y);
}
/* SAM_LISTING_END_0 */
