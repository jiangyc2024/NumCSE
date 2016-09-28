//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>
#include <iomanip>
#include <cmath>


int main() {

    // Lambda function, call with sinh(x)
    /* SAM_LISING_BEGIN_1 */
    auto sinh = [] (double x) {
        double t = std::exp(x);
        return .5 * (t - 1./t);
    };
    /* SAM_LISING_END_1 */

    /* SAM_LISING_BEGIN_2 */
    // TODO: compute relative error for $10^{-k}, k  =1,...10$
    /* SAM_LISING_END_2 */
}
