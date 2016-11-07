#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>

#include <Eigen/Dense>

#include <figure.hpp>

using namespace Eigen;

/*!
 * \brief The LinearInterpolant class
 */
/* SAM_LISTING_BEGIN_0 */
class LinearInterpolant {
public:

    /*!
     * \brief LinearInterpolant builds interpolant from data.
     * Sort the array for the first time:
     * the data is not assumed to be sorted
     * sorting is necessary for binary search
     * \param TODO
     */
    LinearInterpolant(/* TODO: pass data here */);

    /*!
     * \brief operator () Evaluation operator.
     * Return the value of $I$ at $x$, i.e. $I(x)$.
     * Performs bound checks (i.e. if $x < t_0$ or $x >= t_n$ )
     * \param x Value $x \in \mathbf{R}$.
     * \return Value $I(x)$.
     */
    double operator() (double x);
private:
    // TODO: your data there
};
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
LinearInterpolant::LinearInterpolant(/* TODO: pass data here */) {
    // TODO: construct your data here
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double LinearInterpolant::operator() (double x) {
    // TODO: your data there
    return 0.;
}
/* SAM_LISTING_END_2 */

int main(void) {
    // Test the class with the basis with nodes (-1,1,2,4)
    // and interpolant with values (-1,2,3,4).
    LinearInterpolant I;

    int n = 100;
    VectorXd x = VectorXd::LinSpaced(n, -1, 4);
    VectorXd y = x.unaryExpr(
                [&I] (double x) { return I(x); }
    );

    mgl::Figure fig;
    fig.title("Piecewise linear interpolation polynomial.");
    fig.plot(x, y, "b").label("I(x)");
    fig.xlabel("$x$");
    fig.ylabel("$y$");
    fig.legend(0, 1);
    fig.save("plp.eps");
    fig.save("plp.png");
}
