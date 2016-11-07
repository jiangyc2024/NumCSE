#include <iostream>

#include <Eigen/Dense>

#include <figure.hpp>

#include "polyfit.hpp"
#include "polyval.hpp"

using namespace Eigen;

/*!
 * \brief dipoleval
 * \param t
 * \param y
 * \param x
 * \return
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd dipoleval(const VectorXd & t,
                   const VectorXd & y,
                   const VectorXd & x) {
    assert(t.size() = y.size()
           && "t and y must have same size!");

    VectorXd ret;
    ret.resizeLike(x);

#if SOLUTION
    for(int i = 0; i < x.size(); ++i ) {
        VectorXd p(y);
        VectorXd dP = VectorXd::Zero(y.size());

        for(int im = 1; im < y.size(); ++im ) {
            for(int i0 = im-1; i0 >= 0; --i0) {
                dP(i0) = (
                            p(i0+1) + (x(i)-t(i0))*dP(i0+1)
                            - p(i0) - (x(i)-t(im))*dP(i0)
                         ) / (
                            t(im) - t(i0)
                         );
                p(i0) = (
                            (x(i)-t(i0))*p(i0+1)
                          - (x(i)-t(im))*p(i0)
                        ) / (
                            t(im) - t(i0)
                        );
            }
        }

        ret(i) = dP(0);
    }
#else // TEMPLATE
    // TODO Evaluate derivative of interpolating polynomial using AN-scheme
#endif
    return ret;
}
/* SAM_LISTING_END_0 */

/*!
 * \brief dipoleval_alt
 * \param t
 * \param y
 * \param x
 * \return
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd dipoleval_alt(const VectorXd & t,
                       const VectorXd & y,
                       const VectorXd & x) {
    assert(t.size() = y.size()
           && "t and y must have same size!");

    int n = y.size();

#if SOLUTION
    VectorXd P = polyfit(t,y,n-1).head(n-1);

    for(int i = 0; i < P.size(); ++i) {
        P(i) *= n-1-i;
    }

    return polyval(P, x);
#else // TEMPLATE
    VectorXd P;
    // TODO Evaluate derivative of interpolating polynomial using polyfit/polyval
    return P;
#endif
}
/* SAM_LISTING_END_1 */

int main() {

    int n = 10;
    int N = 100;

    VectorXd t = VectorXd::LinSpaced(n,0,3);
    VectorXd y = t.unaryExpr([] (double x) { return std::sin(x);});

    VectorXd x = VectorXd::LinSpaced(N,0,3);

    VectorXd P = polyval(polyfit(t,y,n-1), x);
    VectorXd dp = dipoleval(t,y,x);
    VectorXd dp_alt = dipoleval_alt(t,y,x);

    std::cout << "Error: "
              << (dp-dp_alt).norm()
              << std::endl;

    mgl::Figure fig;
    fig.title("Derivative of interpolation polynomial.");
    fig.plot(x, dp, "r").label("p'(t) (AN)");
    fig.plot(x, dp_alt, "g").label("p'(t) (polyfit)");
    fig.plot(t, y, " b+").label("interpolation pts.");
    fig.plot(x, P, "b").label("p(t)");
    fig.xlabel("$x$");
    fig.ylabel("$y$");
    fig.legend(0, 1);
    fig.save("eval_deriv.eps");
    fig.save("eval_deriv.png");
}
