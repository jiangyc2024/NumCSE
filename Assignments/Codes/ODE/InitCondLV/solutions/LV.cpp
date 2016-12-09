#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "ode45.hpp"

using namespace Eigen;

/*!
 * \brief Compute the maps Phi and W at time T.
 * Use initial data given by u0 and v0.
 * \param[in] u0 First component.
 * \param[in] v0 Second component.
 * \param[in] T Final time.
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<Vector2d, Matrix2d> PhiAndW(double u0,
                                      double v0,
                                      double T) {
    std::pair<Vector2d, Matrix2d> PaW;
    auto f = [] (const VectorXd & w) {
        Eigen::VectorXd temp(6);
        temp(0) = (2. - w(1))*w(0);
        temp(1) = (w(0) - 1.)*w(1);
        temp(2) = (2. - w(1))*w(2) - w(0)*w(3);
        temp(3) = w(1)*w(2) + (w(0) - 1.)*w(3);
        temp(4) = (2. - w(1))*w(4) - w(0)*w(5);
        temp(5) = w(1)*w(4) + (w(0) - 1.)*w(5);
        return temp;
    };

    Eigen::VectorXd w0(6);
    w0 << u0, v0, 1., 0, 0, 1.;

    // Construct ode solver with r.h.s
    ode45<Eigen::VectorXd> O(f);
    // Set options
    O.options.rtol = 1e-14;
    O.options.atol = 1e-12;
    // Solve ODE
    auto sol = O.solve(w0, T);
    // Extract needed component
    VectorXd wT = sol.back().first;

    PaW.first  << wT(0), wT(1);
    PaW.second << wT(2), wT(4), wT(3), wT(5);
    return PaW;
}
/* SAM_LISTING_END_1 */

int main(){
    /* SAM_LISTING_BEGIN_2 */
    // Initial guess
    Vector2d y;
    y << 3, 2;
    // Period we require
    double T = 5;

    // Compute Phi and W from first guess
    std::pair<Vector2d, Matrix2d> PaW = PhiAndW(y(0), y(1), T);
    // Check error
    Vector2d F = PaW.first - y;
    Matrix2d DF;

    // Until we are happy with our approximation
    while (F.norm() > 1e-5) {
        // Test current guess
        PaW = PhiAndW(y(0), y(1), T);
        // Find out error (we want to find a zero of the error)
        F = PaW.first - y;
        // Find out Jacobian
        DF = PaW.second - MatrixXd::Identity(2,2);
        // Use newton iteration
        y = y - DF.lu().solve(F);
    }

    std::cout << "The obtained initial condition is: "
         << std::endl << y << std::endl;
    PaW = PhiAndW(y(0), y(1), 100);

    std::cout << "y(100) = " << std::endl
              << PaW.first << std::endl;
    /* SAM_LISTING_END_2 */
}
