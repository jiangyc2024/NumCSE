//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <limits>
#include <cmath>

double myfunction(double x) {
    double log2=0.693147180559945;
    double y=0;
    while(x>std::sqrt(2)){x/=2; y+=log2;}
    while(x<1./std::sqrt(2)){x*=2; y-=log2;}
    double z=x-1;
    double dz=x*std::exp(-z)-1.;
    while(std::abs(dz/z)>std::numeric_limits<double>::epsilon()) {
        z+=dz;dz=x*std::exp(-z)-1;
    }
    return y+z+dz;
}

double myfunction_modified(double x) {
    double log2 = 0.693147180559945;
    double y = 0;
    while(x > std::sqrt(2)) { x/=2; y+=log2; }
    while(x < 1./std::sqrt(2)) { x*=2; y-=log2; }
    double z=x-1;
    double dz=x*std::exp(-z)-1.;
    // TODO: modify  myfunction to perform a fixed number of iterations
    return 0;
}

int main(int argc, char**argv) {
    // x will contain the first argument to the command line
    double x = 2.;
    if(argc > 1) {
        x = std::stod(argv[1]);
    }
    assert(x > 0 && "x must be > 0!");

    std::cout << "x:                  " << x << std::endl
              << "myfunction:         " << myfunction(x) << std::endl
              << "myfunction_modifed: " << myfunction_modified(x) << std::endl
              << std::endl;
}
