///
/// Minimal runner for (8-4) A derivative-free iterative scheme for finding zeros
///

#include <iostream>

#include "quadraticconvergence.hpp"

int main() {
    // Test Steffensen's method
    testSteffensen(); 

    // Compute errors
    std::cout << "Table for convergence order:" << std::endl;
    orderSteffensen();
    
    return 0;
}
