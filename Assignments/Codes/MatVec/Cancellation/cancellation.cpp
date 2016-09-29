#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#if INTERNAL
#include <figure/figure.hpp>
#endif // INTERNAL

using namespace Eigen;

int main() {
#if SOLUTION
    /* SAM_LISTING_BEGIN_1 */
    // Array of values of $h$
    ArrayXd h = ArrayXd::LinSpaced(21, -20, 0.)
        .unaryExpr([] (double i) {
            return std::pow(10., i);
        });
    // Dummy array where to evaluate the derivative (1.2)
    ArrayXd x = ArrayXd::Constant(h.size(), 1.2);

    // Derivative
    ArrayXd g1 = (sin(x +h) - sin(x)) / h; // naive
    ArrayXd g2 = 2 * cos(x+0.5*h) * sin(0.5 * h) / h; // better
    ArrayXd ex = cos(x); // exact

    //// Print error

    // Table header
    std::cout << std::setw(15) << "h"
              << std::setw(15) << "exact"
              << std::setw(15) << "cancellation"
              << std::setw(15) << "error"
              << std::setw(15) << "improved"
              << std::setw(15) << "error" << std::endl;
    for(unsigned int i = 0; i < h.size(); ++i) {
        // Table entries
        std::cout << std::setw(15) << h(i)
                  << std::setw(15) << ex(i)
                  << std::setw(15) << g1(i)
                  << std::setw(15) << std::abs(g1(i) - ex(i))
                  << std::setw(15) << g2(i)
                  << std::setw(15) << std::abs(g2(i) - ex(i)) << std::endl;
    }
    /* SAM_LISTING_END_1 */
#else // TEMPLATE
    // TODO: Compute approximation of the derivative of sin(x)
    // Print the error of each computation
#endif // TEMPLATE

#if INTERNAL && SOLUTION
    // Plot
    mgl::Figure fig;
    fig.setlog(true, true);
    fig.legend();
    fig.title("Error of approximation of f'(x_0)");
    fig.xlabel("h");
    fig.ylabel("| f'(x_0) - g_i(x_0, h) |");
    fig.plot(h.matrix(), (g1-ex).abs().matrix()).label("g_1");
    fig.plot(h.matrix().tail(16), (g2-ex).abs().matrix().tail(16)).label("g_2");
    fig.plot(h.matrix(), h.matrix(), "h;").label("O(h)");
    fig.save("error_cancellation.eps");
#endif // INTERNAL && SOLUTION
}
