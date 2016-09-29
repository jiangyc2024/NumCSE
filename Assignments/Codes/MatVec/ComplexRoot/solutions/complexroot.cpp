#include <iostream>
#include <complex>
#include <cmath>

//! \brief Compute complex root
//! \param[in] w complex number with non negative imaginary parts
//! \return the square root of w with non negative real and imaginary parts
/* SAM_LISTING_BEGIN_1 */
std::complex<double> myroot( std::complex<double> w ) {
    /* SAM_LISTING_BEGIN_1 */
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
