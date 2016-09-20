#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"

#include "strassen.hpp"

int main() {
    /* SAM_LISTING_BEGIN_1 */
    // Store seed for rng
    unsigned seed = (unsigned) time(0);
     // seed random number generator
    srand(seed);

    // Minimum number of repetitions
    unsigned int repetitions = 10;

    // TODO: time Strassen and Eigen matrix multiplication
    // repeat 10 times and output the minimum runtime
    // use 3 digits and scientific notation
    /* SAM_LISTING_END_1 */
}
