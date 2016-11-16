#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>

#include <Eigen/Dense>

#include <figure.hpp>

using namespace Eigen;

/*!
 * \brief trigIpL
 * \param n
 * \return
 */
/* SAM_LISTING_END_3 */
double trigIpL(std::size_t n) {
    double ret = 0;

#if SOLUTION
    ArrayXd t = ArrayXd::LinSpaced(1e4, 0, 1);

    ArrayXd sint = t.unaryExpr([] (double t) {
        return std::sin(M_PI*t);
    });
    for(unsigned int j = 0; j <= 2*n; ++j) {
        auto bj = [j] (double t) {
            return std::sin(2.*M_PI*(j+.5)*t);
        };

        auto trim_nans = [] (double t) {

            return isnan(t) ? 0 : t;
        };

        ret += (t.unaryExpr(bj) / sint)
                .unaryExpr(trim_nans)
                .cwiseAbs()
                .maxCoeff();

        std::stringstream title, name, legend;
        title << "b_j, j = " << j;
        name << "b_j, j = " << j;
        legend << "b_j, j = " << j;
#if INTENRAL
        mgl::Figure fig;
        fig.title(title.str().c_str());
        fig.xlabel("t");
        fig.ylabel("y");
        fig.plot(t, t.unaryExpr(bj) / sint, "r").label(legend.str().c_str());
        fig.legend();
        fig.save(name.str().c_str());
#endif // INTERNAL
    }
    return ret / 2. / (n + 1/2.);
#else // TEMPLATE

#endif // TEMPLATE
}
/* SAM_LISTING_END_0 */

int main() {
    const int s = 11;

    std::cout << std::setw(s) << "2^k"
              << std::setw(s) << "lambda(k)" << std::endl;

    for(unsigned int i = 1 << 2; i < (1 << 15); i = i << 1) {
        std::cout << std::setw(s) << i
                  << std::setw(s) << trigIpL(i) << std::endl;

    }
}
