#include <iostream>
#include <limits>
#include <cmath>
#include <cassert>

/* SAM_LISTING_BEGIN_1 */
double myfunction(double x) {
    double log2=0.693147180559945;
    double y=0;
    while(x>2*std::sqrt(2)){x/=2; y+=log2;} // \Label[line]{cq:1}
    while(x<std::sqrt(2)){x*=2; y-=log2;} // \Label[line]{cq:2}
    double z=x-1; // \Label[line]{cq:3}
    double dz=x*std::exp(-z)-1.0;
    while(std::abs(dz/z)>std::numeric_limits<double>::epsilon()) {
        z+=dz; dz=x*std::exp(-z)-1.0;
    }
    return y+z+dz; // \Label[line]{cq:4}
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double myfunction_modified(double x) {
    double log2 = 0.693147180559945;
    double y = 0;
    while(x > 2*std::sqrt(2)) { x/=2; y+=log2; }
    while(x < std::sqrt(2)) { x*=2; y-=log2; }
    double z=x-1;
    double dz=x*std::exp(-z)-1.;
    // TODO: modify myfunction to perform a fixed number of iterations
    return 0;
}
/* SAM_LISTING_END_2 */

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
