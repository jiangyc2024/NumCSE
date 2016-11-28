#include <iostream>
#include <limits>
#include <cmath>
#include <cassert>

/*!
 * \brief myfunction Implement log(x) using Newton method
 * \param x A numer > 0
 * \return log(x)
 */
/* SAM_LISTING_BEGIN_1 */
double myfunction(double x) {
    // Step 1: get close to the zero
    double log2=0.693147180559945;
    double y=0;
    while(x>std::sqrt(2)){x/=2; y+=log2;}
    while(x<1./std::sqrt(2)){x*=2; y-=log2;}
    // Step 2: Newton iteration up to terminatio
    double z=x-1;
    double dz=x*std::exp(-z)-1.;
    while(std::abs(dz/z)>std::numeric_limits<double>::epsilon()) {
        z+=dz;dz=x*std::exp(-z)-1;
    }
    return y+z+dz;
}
/* SAM_LISTING_BEGIN_1 */

/*!
 * \brief myfunction_modified Implement log(x) using Newton method
 * Termination criteria choosed a priori, fixed number of iterations.
 * \param x A numer > 0
 * \return log(x)
 */
/* SAM_LISTING_BEGIN_2 */
double myfunction_modified(double x) {
    double log2 = 0.693147180559945;
    double y = 0;
    while(x > std::sqrt(2)) { x/=2; y+=log2; }
    while(x < 1./std::sqrt(2)) { x*=2; y-=log2; }
    double z=x-1;
    double dz=x*std::exp(-z)-1.;
    // Step 2: modifed to perform a fixed number of iterations
    double e0 = z - std::log(x);
    double eps = std::numeric_limits<double>::epsilon();
    double k= (std::log(-std::log(eps))-std::log(-std::log(std::abs(e0))))/std::log(2);
    for (int i=1; i<k; ++i) {
        z = z+dz;
        dz = x * std::exp(-z) - 1;
    }
    return y+z+dz;
}
/* SAM_LISTING_BEGIN_2 */

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
              << "log(x):             " << std::log(x)
              << std::endl;
}
