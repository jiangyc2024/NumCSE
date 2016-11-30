#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>

#if INTERNAL
#include <figure.hpp>
#endif // INTERNAL

using namespace Eigen;

/* @brief Newton's method to approximate $x^{(0)}$
 * @param[in] x0_ Initial guess
 * @param[out] x0 Final estimation of $x^{(0)}$, given convergence of Newton's method
 */
/* SAM_LISTING_BEGIN_0 */
double newton_arctan(double x0_) {

    double x0 = x0_;

#if SOLUTION
    double upd = 1;
    double eps = std::numeric_limits<double>::epsilon();

    while(std::abs(upd) > eps) {

        double x1 = x0 - (2*x0-(1+x0*x0)*std::atan(x0)) / (1-2*x0*std::atan(x0));
//        double x1 = (-x0+(1-x0*x0)*std::atan(x0)) / (1-2*x0*std::atan(x0));

        upd = (x1 - x0) / x1;
        x0 =  x1;
        std::cout << "x0 = " << x1 << ", accuracy = " << upd << std::endl;
    }
#else // TEMPLATE
// TODO: approximate $x^{(0)}$ using Newton's method
#endif // TEMPLATE

    return x0;
}
/* SAM_LISTING_END_0 */

int main() {
	// Initialization
    double x0_ = 2; // Initial guess

    // Netwon's method
    double x0 = newton_arctan(x0_);
    double x1 = x0 - std::atan(x0)*(1+x0*x0);
    double x2 = x1 - std::atan(x1)*(1+x1*x1);

#if INTERNAL
    // Figure
    VectorXd X = VectorXd::LinSpaced(401,-2,+2);
    VectorXd atan(X.size());
    for(unsigned i=0; i<X.size(); ++i) {
        atan(i) = std::atan(X(i));
    }
    VectorXd g = X.array().square();
    g = 2*X - (VectorXd::Ones(X.size()) + g).cwiseProduct(atan);
    VectorXd X_1(5); X_1 << x0, x1, x1, x2, x2;
    VectorXd Y_1(5); Y_1 << std::atan(x0), 0, std::atan(x1), 0, std::atan(x2);
    VectorXd X_2(2); X_2 << x0, x1;
    VectorXd Y_2(2); Y_2 << 0, 0;

    mgl::Figure fig;
    fig.title("Newton critical iteration");
    fig.ranges(-2, +2, -1.2, +1.2);
    fig.plot(X, atan, "k-").label("arctan");
    fig.plot(X, g, "r=").label("g");
    fig.plot(X_1, Y_1, "b-").label("Newton critical iteration");
    fig.plot(X_2, Y_2, "o ");
    fig.fplot("0", "k-");
    fig.legend(0, 1);
    fig.save("NewtonArctan");
#endif // INTERNAL
}
