#include <iostream>
#include <iomanip>
#include <cmath>

#if SOLUTION
#include "sinh_stable.hpp"
#endif // SOLUTION

int main() {

    // Lambda function, call with sinh(x)
    /* SAM_LISING_BEGIN_1 */
    auto sinh = [] (double x) {
        double t = std::exp(x);
        return .5 * (t - 1./t);
    };
    /* SAM_LISING_END_1 */

    /* SAM_LISING_BEGIN_2 */
#if SOLUTION
    std::cout << "--> Realitve error of varoious implementations:" << std::endl;
    std::cout << std::setw(10) << "k"
              << std::setw(15) << "own"
              << std::setw(15) << "std"
              << std::setw(15) << "tylor"
              << std::setw(15) << "err"
              << std::setw(15) << "err (taylor)"
              << std::endl;
    for( int k = 1; k <= 10; ++k ) {
        // Value x
        double x = std::pow<double>(10., -k);

        // Our own lambda function
        double mySinh = sinh(x);
        // The "standard" sinh
        double stdSinh = std::sinh(x);
        // The stable sinh (look at the file "sinh_stable.hpp" if
        // interested in advanced C++)
        double taylorSinh = taylor_sinh<3>(x);

        // Relative error
        double err = std::abs(mySinh - stdSinh) / std::abs(stdSinh);
        double err_t = std::abs(taylorSinh - stdSinh) / std::abs(stdSinh);

        // Print error
        std::cout << std::setw(10) << k
                  << std::setw(15) << mySinh
                  << std::setw(15) << stdSinh
                  << std::setw(15) << taylorSinh
                  << std::setw(15) << err
                  << std::setw(15) << err_t << std::endl;
    }

    // We compute and print the error bound computed using the
    // remainder of the Taylor expansion
    double x = 1e-3;
    std::cout << "--> Relative error bound for x = " << x << ":"
              << std::endl;
    std::cout << std::setw(10) << "n"
              << std::setw(15) << "error bound" << std::endl;
    std::cout << std::setw(10) << "1"
              << std::setw(15) << error_bound<1>(x) << std::endl;
    std::cout << std::setw(10) << "2"
              << std::setw(15) << error_bound<2>(x) << std::endl;
    std::cout << std::setw(10) << "3"
              << std::setw(15) << error_bound<3>(x) << std::endl;
    std::cout << std::setw(10) << "4"
              << std::setw(15) << error_bound<4>(x) << std::endl;

#else // TEMPLATE
    // TODO: compute relative error for $10^{-k}, k  =1,...10$
#endif // TEMPLATE
    /* SAM_LISING_END_2 */
}
