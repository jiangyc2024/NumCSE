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
double newton_arctan(double x0) {

    double x0_ = x0;

#if SOLUTION
    double upd = 1;
    double eps = std::numeric_limits<double>::epsilon();

    while(std::abs(upd) > eps) {

        double x1 = x0_ - (2*x0_-(1+x0_*x0_)*std::atan(x0_)) / (1-2*x0_*std::atan(x0_));
//        double x1 = (-x0_+(1-x0_*x0_)*std::atan(x0_)) / (1-2*x0_*std::atan(x0_));

        upd = (x1 - x0_) / x1;
        x0_ =  x1;
        std::cout << "x0 = " << x1 << ", accuracy = " << upd << std::endl;
    }
#else // TEMPLATE
// TODO:
#endif // TEMPLATE

    return x0_;
}
/* SAM_LISTING_END_0 */

int main() {
	// Initialization
    double x0  = 2; // Initial guess

    // Netwon's method
    double x0_ = newton_arctan(x0);

//    figure;
//    x1 = x0-atan(x0)*(1+x0^2);   x2 = x1-atan(x1)*(1+x1^2);
//    X=[-2:0.01:2];
//    plot(X, atan(X),'k',...
//        X,  2*(X)-(1+(X).^2).*atan((X)),'r--',...
//        [x0, x1, x1, x2, x2], [atan(x0), 0, atan(x1), 0, atan(x2)],...
//        [x0,x1],[0,0],'ro',[-2,2], [0,0],'k','linewidth',2);
//    legend('arctan', 'g', 'Newton critical iteration');axis equal;
//    print -depsc2 'ex_NewtonArctan.eps'
}
