//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>
#include <complex>
#include <cmath>

//! \brief Compute complex root
//! \param[in] w complex number with non negative imaginary parts
//! \return the square root of w with non negative real and imaginary parts
std::complex<double> myroot( std::complex<double> w ) {
    double x,y;
    double u = w.real();
    double v = w.imag();

    if (v==0) return std::sqrt(u);

    if (u > 0) {
        x = std::sqrt((std::sqrt(u*u+v*v)+u)/2.);
        y = v/(2*x);
    } else {
        y = std::sqrt((std::sqrt(u*u+v*v)-u)/2.);
        x = v/(2*y);
    }

    return std::complex<double> (x,y);
}

// Test the implementation
int main() {
    std::cout << "--> Correctness test:" << std::endl;

    std::complex<double> w(1e20, 5);
    std::cout << "The square root of " << w << " is " << myroot(w) << std::endl ;
    std::cout << "The correct square root of " << w << " is " << std::sqrt(w) <<std::endl <<std::endl;

    w = std::complex<double>(-5, 1e20);
    std::cout << "The square root of " << w << " is " << myroot(w) << std::endl ;
    std::cout << "The correct square root of " << w << " is " << std::sqrt(w) <<std::endl <<std::endl;
}
