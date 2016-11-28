#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>

using namespace Eigen;

/* @brief
 * @param[in]
 * @param[out]
 */
/* SAM_LISTING_BEGIN_0 */
template <class function>
VectorXd steffensen(function&& f, double x0) {

    VectorXd x(1); x(0) = x0;

#if SOLUTION
    double upd = 1;
    double eps = std::numeric_limits<double>::epsilon();

    while(std::abs(upd) > eps) {

        fx = f(x.tail(1)); // Only 2 evaluations of $f$ at each step
        if(fx != 0) {
            upd = fx*fx / (f(x.tail(1)+fx)-fx);
            x_new = x.tail(1) - upd;
            x.conservativeResize(x.size()+1);
            x.tail(1) = x_new;
        } else {
            upd = 0;
        }
    }
#else // TEMPLATE
// TODO:
#endif // TEMPLATE

    return x;
}
/* SAM_LISTING_END_0 */

int main() {
	// Initialization
    using Vector = Eigen::VectorXd;
    // Different definitions of $f$ with the same zero:
    auto f = [] (double x) { return x*std::exp(x)-1; };
//    auto f = [] (double x) { return std::exp(x)-1/x; };
//    auto f = [] (double x) { return x-std::exp(-x); };
    double x0 = 1;
    double x_star = 0.567143290409784; // From Matlab: x_star = fzero(f,x0);

    // Steffensen's method
    VectorXd x = steffensen(f, x0);

    // Compute errors
    double residual = f(x);
    double err = std::abs(x-x_star);
    double log_err = std::log(err);
    VectorXd ratios(log_err.size()-2);
    for(unsigned i=0; i<ratios.size(); ++i) {
        ratios(i) = (log_err(i+2) - log_err(i+1)) /
                    (log_err(i+1) - log_err(i));
    }

    // Print output
    std::cout << "x" << "\t" << "err" << "\t" << "residuals" <<" \t" << "ratios" << std::endl;
    std::cout << x(0) << "\t" << err(0) << "\t" << residuals(0) << "\t" << ratios(0) << std::endl;
    std::cout << x(1) << "\t" << err(1) << "\t" << residuals(1) << "\t" << ratios(1) << std::endl;
    for(unsigned i=0; i<ratios.size(); ++i) {
        std::cout << x(i+2) << "\t" << err(i+2) << "\t" << residuals(i+2) << "\t" << ratios(i) << std::endl;
    }
}
