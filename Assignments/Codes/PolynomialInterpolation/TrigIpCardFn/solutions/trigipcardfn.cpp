#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>

#include <Eigen/Dense>

#include <figure.hpp>

#include "trigpolyvalequid.hpp"

using namespace Eigen;

/*!
 * \brief plot_basis Plot the shifted basis polynomials.
 * \param n $2*n+1$ will be the number of basis polynomials.
 */
void plot_basis(int n) {
    // Mesh size
    const int M = 1e3;

    // Basis vector e_1
    ArrayXd e = ArrayXd::Zero(2*n+1);
    e(0) = 1;
    VectorXd y;
    trigpolyvalequid(e, 1e3, y);

    ArrayXd t = ArrayXd::LinSpaced(M, 0, 1);

    // Shift function right a bit
    ArrayXd y_shift(M);
    unsigned int h = M / (2*n+1);
    y_shift << y.tail(h), y.head(M - h);

    mgl::Figure fig;
    fig.title("b_0(t)");
    fig.xlabel("t");
    fig.ylabel("y");
    fig.plot(t, y_shift, "r").label("b_0(t)");
    fig.legend();
    fig.save("b0_n");
}

/*!
 * \brief trigIpL Compute $\lambda(n)$.
 *
 * \param[in] n $2*n+1$ will be the number of basis polynomials.
 * \return Value $\lambda(n)$.
 */
/* SAM_LISTING_BEGIN_0 */
double trigIpL(std::size_t n) {
    double ret = 0;

    ArrayXd t = ArrayXd::LinSpaced(1e4, 0, 1);

    ArrayXd sint = t.unaryExpr([] (double t) {
        return std::sin(M_PI*t);
    });
    for(unsigned int j = 0; j <= 2*n; ++j) {
        auto bj = [j] (double t) {
            return std::sin(2.*M_PI*(j+.5)*t);
        };

        auto trim_nans = [] (double t) {

            return std::isnan(t) ? 0 : t;
        };

        ret += (t.unaryExpr(bj) / sint)
                .unaryExpr(trim_nans)
                .cwiseAbs()
                .maxCoeff();

    }
    return ret / 2. / (n + 1/2.);
}
/* SAM_LISTING_END_0 */

int main() {

    /// PART 1
    const int n = 5;
    plot_basis(n);

    /// PART 2
    const int s = 11;

    std::cout << std::setw(s) << "2^k"
              << std::setw(s) << "lambda(k)" << std::endl;

    for(unsigned int i = 1 << 2; i < (1 << 15); i = i << 1) {
        double l = trigIpL(i);
        std::cout << std::setprecision(3)
                  << std::setw(s) << i
                  << std::setw(s) << l << std::endl;

    }
}
