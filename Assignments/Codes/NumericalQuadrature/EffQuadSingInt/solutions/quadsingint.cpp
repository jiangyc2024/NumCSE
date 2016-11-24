#include <iomanip>
#include <iostream>
#include <cmath>

// Include function "gauleg" for the computation of weights/nodes of Gauss quadrature
#include "gauleg.hpp"
#include "polyfit.hpp"
#include "polyval.hpp"

// Choose between substitution 1. s = sqrt(1 \pm t), 2. t = cos(s)
#define METHOD 1
// Set to true to enable test  with monomials
#define MONOMIAL_TEST   false

//! @brief Compute integral $\int_{-1}^1 \sqrt(1-t^2) f(t) dt$ using tranfromation
//! @param[in] n number of Gauss nodes (evaluate f at 2*n points)
//! @return value of integral
/* SAM_LISTING_BEGIN_1 */
template <class Function>
double quadsingint(const Function& f, const unsigned n) {
    double I = 0.;
    // Method selects one of the two possible solutions
#if METHOD == 1 // $s = \sqrt(1 \pm t)$
    // Will use same node twice
    QuadRule Q = gauleg(n);

    for(unsigned i = 0; i < n; ++i) {
        // Transform nodes
        double x = (Q.nodes(i) + 1.) / 2.;
        // Weights
        double w = Q.weights(i) * x * x * std::sqrt(2. - x * x);
        // Symmetric summation
        I += w * (f(x*x - 1) + f(-x*x + 1));
    }
#elif METHOD == 2 // $t = cos(s)$
    // We are actually using twice the number of nodes than Method 1
    QuadRule Q = gauleg(2*n);
    
    for(unsigned i = 0; i < 2*n; ++i) {
        // Evualuate transformation
        double x = sin(Q.nodes(i) * M_PI_2);
        // Weights
        double w = Q.weights(i) * cos(Q.nodes(i) * M_PI_2) * cos(Q.nodes(i) * M_PI_2);
        I += w * f(x) * M_PI_2;
    }
#endif
    return I;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
int main(int argc, char ** argv) {
    
    // Max num of Gauss points to use (in each direction)
    const unsigned max_N = 25;
    // "Exact" integral
    double I_ex = 0.483296828976607;

    // If you want to test monomials as f, use this data and set MONOMIAL\_TEST to true
#if MONOMIAL_TEST
    std::vector<double> ex = { M_PI_2, 0., M_PI_2 / 4, 0, M_PI_2 / 8, 0, M_PI_2 * 5 / 64};
    assert(argc > 1);
    int n = std::atoi(argv[1]);
    std::cout << "TEST: Monomial of degree:" << n << std::endl;
    auto f = [&n] (double t) { return std::pow(t, n); };
    I_ex = ex[n];
#else // END MONOMIAL\_TEST
    // Test function
    auto f = [] (double t) { return 1. / (2. + std::exp(3*t) ); };
#endif

    std::cout << std::setw(3) << "N"
              << std::setw(15) << "I_approx"
              << std::setw(15) << "error"
              << std::endl;
    double errold, errnew;
    // Observed: exponential convergence (as expected)
    for(unsigned int N = 1; N < max_N; ++N) {
        double I_approx = quadsingint(f, N);

        errold = errnew;
        errnew = std::abs(I_ex - I_approx);

        // Interpolating points
        std::cout << std::setw(3) << N
                     // Value of integral
                  << std::setw(15) << I_approx
                     // Error
                  << std::setw(15) << errnew
                  << std::endl;
    }
    return 0;
}
/* SAM_LISTING_END_2 */
