#include <iostream>

#include "complexroot.hpp"



// Test the implementation
int main() {
    std::cout << "--> Test:" << std::endl;

    std::complex<double> w(1e20, 5);
    std::cout << "The square root of " << w << " is " << myroot(w) << std::endl ;

    w = std::complex<double>(-5, 1e20);
    std::cout << "The square root of " << w << " is " << myroot(w) << std::endl ;
    
    w = std::complex<double>(1e-8, 0);
    std::cout << "The square root of " << w << " is " << myroot(w) << std::endl ;
}
