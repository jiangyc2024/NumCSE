#include <iostream>
#include <iomanip>

#include <Eigen/Dense>


using namespace Eigen;

int main() {
    /* SAM_LISTING_BEGIN_1 */
    ArrayXd h = ArrayXd::LinSpaced(21, -20, 0.)
        .unaryExpr([] (double i) {
            return std::pow(10., i);
        });
    ArrayXd x = ArrayXd::Constant(h.size(), 1.2);

    // Derivative
    ArrayXd g1 = (sin(x +h) - sin(x)) / h; // naive
    ArrayXd g2 = 2 * cos(x+0.5*h) * sin(0.5 * h) / h; // better
    ArrayXd ex = cos(x);

    std::cout << std::setw(15) << "h"
              << std::setw(15) << "exact"
              << std::setw(15) << "cancellation"
              << std::setw(15) << "error"
              << std::setw(15) << "improved"
              << std::setw(15) << "error" << std::endl;
    for(unsigned int i = 0; i < h.size(); ++i) {
        std::cout << std::setw(15) << h(i)
                  << std::setw(15) << ex(i)
                  << std::setw(15) << g1(i)
                  << std::setw(15) << std::abs(g1(i) - ex(i))
                  << std::setw(15) << g2(i)
                  << std::setw(15) << std::abs(g2(i) - ex(i)) << std::endl;
    }
    /* SAM_LISTING_END_1 */

}
