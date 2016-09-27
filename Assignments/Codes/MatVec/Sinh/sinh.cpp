#include <iostream>
#include <iomanip>
#include <cmath>

#include "sinh_stable.hpp"

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
    std::cout << std::setw(10) << "k"
              << std::setw(15) << "own"
              << std::setw(15) << "std"
              << std::setw(15) << "tylor"
              << std::setw(15) << "err"
              << std::setw(15) << "err (taylor)"  << std::endl;
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
#else // TEMPLATE
    // TODO: compute relative error for $10^{-k}, k  =1,...10$
#endif // TEMPLATE
    /* SAM_LISING_END_2 */
}
