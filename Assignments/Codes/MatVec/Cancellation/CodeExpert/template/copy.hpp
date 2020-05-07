#ifndef COPY_HPP
#define COPY_HPP
#include "../solution/cancellation.hpp"

#include <Eigen/Dense>
#include "plot.hpp"

using namespace Eigen;

/* SAM_LISTING_BEGIN_1 */
void sinederv_TEST() {
    // Array of values of $h$
    ArrayXd h = ArrayXd::LinSpaced(21, -20, 0.)
        .unaryExpr([] (double i) {
            return std::pow(10., i);
        });
    // Dummy array where to evaluate the derivative (1.2)
    ArrayXd x = ArrayXd::Constant(h.size(), 1.2);
    
    // Derivatives
    ArrayXd g1, g2, ex;
    
    // TO DO: (2-8.b) For each entry in h, fill in the arrays such that
    // g1 approximates f'(1.2) using the difference quotient directly,
    // g2 approximates f'(1.2) while avoiding cancellation, and
    // ex is the "exact" value of f'(1.2).
    // Hint: You can use Eigen's sin() and cos(), which calculate
    // coefficient-wise results for Array inputs.
    // START
    
    // END
    
    // Print error
    // Table header
    std::cout << std::setw(15) << "h"
              << std::setw(15) << "exact"
              << std::setw(15) << "cancellation"
              << std::setw(15) << "error"
              << std::setw(15) << "improved"
              << std::setw(15) << "error" << std::endl;
    for(unsigned int i = 0; i < h.size(); ++i) {
        std::cout << std::setprecision(6);
        // Table entries
        // TO DO: (2-8.b) Print the i-th row of the table.
        // Keep the order h, ex, g1, |g1 - ex|, g2, |g2 - ex|.
        // Use std::setw(15) as in the table header.
        // You can use std::abs() to calculate absolute values.
        // START
        
        // END
    }
    
    // Plots the contents of the table on a log-log scale.
    py_plot( h, ex, g1, g2 );
}
/* SAM_LISTING_END_1 */

#endif