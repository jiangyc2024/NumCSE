///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <Eigen/Dense>

#include "gramschmidt.hpp"

using namespace Eigen;

/* SAM_LISTING_BEGIN_1 */
int main () { /* SAM_SOLUTION_BEGIN */
    // Ortho test
    unsigned int n = 9;
    MatrixXd A = MatrixXd::Random(n,n);
    MatrixXd Q = gramschmidt( A );

    // Output should be idenity matrix
    std::cout << Q*Q.transpose() << std::endl;
    return 0; /* SAM_SOLUTION_END */
}
/* SAM_LISTING_END_1 */