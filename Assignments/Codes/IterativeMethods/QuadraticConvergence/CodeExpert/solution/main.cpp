#include <cmath>
#include <iostream>
#include <limits>


#include "quadraticconvergence.hpp"
int main() {

    // test Steffensen's method
    testSteffensen(); 
    
    // test Logger class
    Logger<double> test_logger; 
    
    test_logger(2.7);
    test_logger(2.0);
    test_logger(2.5);
    std::cout << "The Logger test_logger contains the values:" << std::endl; 
    test_logger.print_log();
    std::cout << std::endl;
    
    // Compute errors
    std::cout << "Table for convergence order:" << std::endl;
    orderSteffensen();
    
    return 0;
}
