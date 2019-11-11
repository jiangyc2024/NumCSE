#include "quadinf.hpp"

int main() {
    
    int n = 10;
    double I = quadinf(n, [] (double t) { return 1 / (1 + std::pow(t,2)); });
    
    std::cout << I << std::endl; //Note: exact integral = pi
    
    //plot error
    cvgQuadInf();
    
    return 0;
}
