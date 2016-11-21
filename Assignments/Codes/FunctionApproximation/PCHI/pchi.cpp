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
 * \brief Implements a piecewise cubic Hermite interpolation.
 * Uses equidistant meshes and various methods of slope reconstruction.
 *
 */
class PCHI {
public:
    /*!
     *! \brief Construct the slopes from the data.
     *! Use finite-differences or setting $s'(x_j) = 0$.
     *! \param[in] t Vector of nodes (assumed equidistant and sorted).
     *! \param[in] y Vector of values at nodes $t$.
     *! \param[in] s Flag to set if you want to reconstruct or set the slopes to zero.
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
    assert( n >= 3 && "Need at least three nodes." );
    h = t(1) - t(0);

    //// Reconstruction of the slope,
    switch(s) {
    /* SAM_LISTING_BEGIN_1 */
        /// CASE: assuming $s'(x_j) = 0$ (error: $O(1)$)
        case Slope::Zero:
#if SOLUTION
        c = VectorXd::Zero(n);
#else // TEMPLATE
    // TODO: set the slopes assuming $s'(x_j) = 0$
#endif // TEMPLATE
            break;
    /* SAM_LISTING_END_1 */
    /* SAM_LISTING_BEGIN_2 */
        /// CASE: second order finite differences (error: $O(h^2)$)
        case Slope::Reconstructed:
        default:

#if SOLUTION
            // First order alternative:
//                c(0) = ( y(1) - y(0) ) / h;
            c(0) = ( -1*y(2) + 4*y(1) - 3*y(0) ) / 2 / h;
            for(int i = 1; i < n-1; ++i) {
                c(i) = ( y(i+1) - y(i-1) ) / 2 / h;
            }

            // First order alternative:
//                c(n-1) = ( y(n-1) - y(n-2) ) / h;
            c(n-1) = ( 3*y(n-1) - 4*y(n-2) + 1*y(n-3) ) / 2 / h;

#else // TEMPLATE
        // TODO: set the slopes using finite difference reconstruction
#endif // TEMPLATE
            break;
    /* SAM_LISTING_END_2 */
    }
}

/* SAM_LISTING_BEGIN_3 */
VectorXd PCHI::operator() (const VectorXd & x) const {
    VectorXd ret(x.size());

#if SOLUTION
    int i_star = 0;
    double tmp, t1, y1, y2, c1, c2, a1, a2, a3;
    for(int j = 0; j < x.size();) {
        // Stores the current interval index and some temporary variable
        // Find the interval for $x_j$ (porting of 'hermloceval' Matlab code 3.4.6)
        // while assuming that $x$ values are sorted and within $[t_{\min},t_{\max}]$
        if( t(i_star) <= x(j) && x(j) <= t(i_star+1)) {
            t1 = t(i_star);
            y1 = y(i_star);
            y2 = y(i_star+1);
            c1 = c(i_star);
            c2 = c(i_star+1);
            a1 = y2 - y1;
            a2 = a1 - h*c1;
            a3 = h*c2 - a1 - a2;

            // Compute $s(x(j))$
            tmp = ( x(j) - t1 ) / h;
            ret(j) = y1 + ( a1 + ( a2+a3*tmp ) * ( tmp - 1. ) ) * tmp;

            ++j;
        } else {
            // Otherwise, go to next interval
            ++i_star;
            // Terminate if outside any interval
            if(i_star >= t.size() - 1) break;
            continue;
        }
    }
#else // TEMPLATE
    // TODO: implement evaluation operator, returning $s(x)$.
#endif // TEMPLATE

    return ret;
}
/* SAM_LISTING_END_3 */

int main() {
    // Interpoland
    auto f = [] (double x) { return 1. / (1. + x*x); };
    // auto f = [] (double x) {return cos(x); };

    double a = 5; // Interval bounds will be (-a,a)
    int M = 1000; // Number of  points in which to evaluate the interpoland

    /* SAM_LISTING_BEGIN_4 */
    // Precompute values at which evaluate f
    VectorXd x = VectorXd::LinSpaced(M, -a, a);
    VectorXd fx = x.unaryExpr(f);

    // Store error and number of nodes
    std::vector<double> N, err_reconstr, err_zero;

    for(int i = 4; i <= 512; i = i << 1) {
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

        std::cout << i << " " <<  err_zero.back() << " " <<  err_reconstr.back() << std::endl;

        // See how interpolant looks
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
    /* SAM_LISTING_END_4 */
}
