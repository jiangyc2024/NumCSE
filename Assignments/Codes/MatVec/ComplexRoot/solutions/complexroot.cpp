#include <iostream>
#include <complex>
#include <cmath>

//! \brief Compute complex root
//! \param[in] w complex number with non negative imaginary parts
//! \return the square root of w with non negative real and imaginary parts
/* SAM_LISTING_BEGIN_1 */
std::complex<double> myroot( std::complex<double> w ) {
    double x,y;
    double u = w.real();
    double v = w.imag();

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

    return std::complex<double> (x,y);
}
/* SAM_LISTING_END_1 */

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
