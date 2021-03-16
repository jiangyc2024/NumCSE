#include <cmath>
#include <iostream>

#include "newtarctan.hpp"

int main() {
  
    double x0 = newton_arctan();
    // Netwon's method on arctan
    double x1 = x0 - std::atan(x0)*(1+x0*x0);
    double x2 = x1 - std::atan(x1)*(1+x1*x1);
    
    std::cout << "x0 = " << x0 << std::endl;
    std::cout << "x1 = " << x1 << std::endl;
    std::cout << "x2 = " << x2 << std::endl;
}

