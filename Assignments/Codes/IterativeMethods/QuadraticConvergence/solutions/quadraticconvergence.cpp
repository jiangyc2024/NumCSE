#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>


using namespace Eigen;

/* @brief Steffensen's method
 * @param[in] f Function handle
 * @param[in] x0 Initial guess
 * @param[out] x All estimations returned by Steffensen's iterations until convergence
 */
/* SAM_LISTING_BEGIN_0 */
template <class function>
VectorXd steffensen(function&& f, double x0) {

    VectorXd x(1); x(0) = x0;

    double upd = 1;
    double eps = std::numeric_limits<double>::epsilon();

    while(std::abs(upd) > eps) {

        double fx = f(x(x.size()-1)); // Only 2 evaluations of $f$ at each step
        if(fx != 0) {
            upd =  fx*fx / (f(x(x.size()-1)+fx)-fx);
            double x_new = x(x.size()-1) - upd;
            x.conservativeResize(x.size()+1);
            x(x.size()-1) = x_new;
        } else {
            upd = 0;
        }
    }

    return x;
}
/* SAM_LISTING_END_0 */

int main() {
    // Initialization
    // Different definitions of $f$ with the same zero:
    auto f = [] (double x) { return x*std::exp(x)-1; };
//    auto f = [] (double x) { return std::exp(x)-1/x; };
//    auto f = [] (double x) { return x-std::exp(-x); };
    double x0 = 1;
    double x_star = 0.567143290409784; // From Matlab: x_star = fzero(f,x0);

    // Steffensen's method
    VectorXd x = steffensen(f, x0);

    // Compute errors
    unsigned n = x.size();
    VectorXd residuals(n), errs(n), log_errs(n);
    for(unsigned i=0; i<n; ++i) {
        residuals(i) = f(x(i));
        errs(i) = std::abs(x(i)-x_star);
        log_errs(i) = std::log(errs(i));
    }
    VectorXd ratios = VectorXd::Zero(n);
    for(unsigned i=2; i<n; ++i) {
        ratios(i) = (log_errs(i)   - log_errs(i-1)) /
                    (log_errs(i-1) - log_errs(i-2));
    }

    // Print output
    std::cout << "x" << "\t" << "errors" <<" \t" << "ratios" << std::endl;
    for(unsigned i=0; i<n; ++i) {
        std::cout << x(i) << "\t" << errs(i) << "\t" << ratios(i) << std::endl;
    }

}
