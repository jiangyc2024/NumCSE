#include <iostream>
#include <complex>
#include <cmath>

//! \brief Compute complex root
//! \param[in] w complex number with non negative imaginary parts
//! \return the square root of w with non negative real and imaginary parts
/* SAM_LISTING_BEGIN_1 */
std::complex<double> myroot( std::complex<double> w ) {
    /* SAM_LISTING_BEGIN_1 */
    // TODO: Implement complex root without cancellation risk
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
