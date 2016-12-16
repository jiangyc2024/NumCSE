#include "rkintegrator.hpp"

//! \file stabrk.cpp Solve prey/predator model with RK-SSM method

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector2d predprey(Eigen::Vector2d y0, double T, unsigned N)
{
    double h = T / N;
    Eigen::Vector2d y_ = y0;

    // TODO: solve the predator/prey model using the 3-stage RK-SSM of the exercise sheet

    return y_;
}
/* SAM_LISTING_END_1 */

int main() {
/* SAM_LISTING_BEGIN_0 */
    // TODO: solve the predator/prey model using class "RKIntegrator"
/* SAM_LISTING_END_0 */
}
