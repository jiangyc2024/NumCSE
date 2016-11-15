#include <vector>
#include <cassert>

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <figure.hpp>

using namespace Eigen;

// Flag for slope reconstruction type
enum class Slope { Zero, Reconstructed };

/* SAM_LISTING_BEGIN_0 */
/*!
 * \brief Implements a piecewise cubic Herite interpolation.
 * Uses equidistant meshes and various methods of slope reconstruction.
 *
 */
class PCHI {
public:
    /*!
     *! \brief Construct the slopes frome the data
     *! Either using dinite-differences or assuming $f'(x_j) = 0$
     *! \param[in] t vector of nodes (assumed equidistant and sorted)
     *! \param[in] y vector of values at nodes $t$
     *! \param[in] s Flag to set if you want to reconstruct or set slopes to zero
     */
    PCHI(const VectorXd & t,
         const VectorXd & y,
         Slope s = Slope::Reconstructed);

    /*!
     *! \brief Evaluate the intepolant at the nodes $x$.
     *! The input is assumed sorted, unique and inside the interval
     *! \param[in] x The vector of points $x_i$ where to compute $s(x_i)$.
     *! \return Values of interpolant at $x$ (vector)
     */
    VectorXd operator() (const VectorXd & x) const;

private:
    // Provided nodes and values $(t,y)$ to compute spline,
    // Eigen vectors, $c$ contains slopes
    // All have the same size $n$
    VectorXd t, y, c;
    // Difference $t(i)-t(i-1)$
    double h;
    // Size of $t$, $y$ and $c$.
    int n;
};
/* SAM_LISTING_END_0 */

PCHI::PCHI(const VectorXd & t,
           const VectorXd & y,
           Slope s)
    : t(t), y(y), c(t.size()) {
    // Sanity check
    n = t.size();
    assert( n == y.size() && "t and y must have same dimension." );
    assert( n >= 3 && "Need at least two nodes." );
    h = t(1) - t(0);

    //// Reconstruction of the slope,
    switch(s) {
    /* SAM_LISTING_BEGIN_1 */
        /// CASE: assuming $s'(x_j) = 0$ (error: $O(1)$)
        case Slope::Zero:
    // TODO: set the slopes assuming $s'(x_j) = 0$
            break;
    /* SAM_LISTING_END_1 */
    /* SAM_LISTING_BEGIN_2 */
        /// CASE: second order finite differences (error: $O(h^2)$)
        case Slope::Reconstructed:
        default:

        // TODO: set the slopes using finite difference reconstruction
            break;
    /* SAM_LISTING_END_2 */
    }

}

/* SAM_LISTING_BEGIN_3 */
VectorXd PCHI::operator() (const VectorXd & x) const {

    VectorXd ret(x.size());

    // TODO: implement evaluation operator, returning $s(x)$.

    return ret;
}
/* SAM_LISTING_END_3 */

int main() {
    // Interpoland
    auto f = [] (double x) { return 1. / (1. + x*x); };
    // auto f = [] (double x) {return cos(x); };

    double a = 5; // Interval bounds will be (-a,a)
    int M = 1000; // Number of  points in which to evaluate the interpoland

    // Precompute values at which evaluate f
    VectorXd x = VectorXd::LinSpaced(M, -a, a);
    VectorXd fx = x.unaryExpr(f);

    // Store error and number of nodes
    std::vector<double> N, err_reconstr, err_zero;

    for(int i = 2; i <= 512; i = i << 1) {
        // Define subintervals and evaluate f there (find pairs (t,y))
        VectorXd t = VectorXd::LinSpaced(i, -a, a);
        VectorXd y = t.unaryExpr(f);

        // Construct PCHI with zero and reconstructed slopes
        PCHI s_reconstr(t,y),
             s_zero(t, y, Slope::Zero);

        // Evaluate interpolant
        VectorXd s_zero_x = s_zero(x);
        VectorXd s_reconstr_x = s_reconstr(x);

        // Compute infinity norm of error
        err_reconstr.push_back(
                    (s_reconstr_x - fx).lpNorm<Infinity>()
                    );
        err_zero.push_back(
                    (s_zero_x - fx).lpNorm<Infinity>()
                    );
        N.push_back(1. / i);

        // Se how interpolant looks
        if( i == 16 ) {
            {
                mgl::Figure fig;
                fig.title("Interpolant with zero slope");
                fig.xlabel("t");
                fig.ylabel("y");
                fig.plot(x, s_zero_x, "r").label("s_{zero}");
                fig.plot(x, fx, "b--").label("f");
                fig.legend();
                fig.save("p_zero");
            }
            {
                mgl::Figure fig;
                fig.title("Interpolant with reconstructed slope");
                fig.xlabel("t");
                fig.ylabel("y");
                fig.plot(x, s_reconstr_x, "r").label("s_{reconstr}");
                fig.plot(x, fx, "b--").label("f");
                fig.legend();
                fig.save("p_reconstr");
            }
        }
    }

    // Error plot
    mgl::Figure fig;
    fig.title("Error VS no. of nodes");
    fig.setlog(true, true);
    fig.xlabel("h");
    fig.ylabel("max |f(t) - s(t)|");
    fig.fplot("x", "k|").label("O(h)");
    fig.fplot("x^3", "k|").label("O(h^{3})");
    fig.plot(N, err_reconstr, "b").label("s_{reconstr}");
    fig.plot(N, err_zero, "r").label("s_{zero}");
    fig.legend();
    fig.save("pchi_conv");
}
