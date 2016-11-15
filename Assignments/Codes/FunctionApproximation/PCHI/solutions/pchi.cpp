#include <vector>
#include <cassert>

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <figure.hpp>

// Flag for slope reconstruction type
enum class Slope { Zero, Reconstructed };

/* SAM_LISTING_BEGIN_0 */
//! \breif Implements a piecewise cubic Herite interpolation on equidistant meshes.
class PCHI {
public:
    //! \brief Construct the slopes frome the data
    //! Either using dinite-differences or assuming $f'(x_j) = 0$
    //! \param[in] t vector of nodes (assumed equidistant and sorted)
    //! \param[in] y vector of values at nodes $t$
    //! \param[in] s Flag to set if you want to reconstruct or set slopes to zero
    PCHI(const Eigen::VectorXd & t,
         const Eigen::VectorXd & y,
         Slope s = Slope::Reconstructed);

    //! \brief Evaluate the intepolant at the nodes $x$.
    //! The input is assumed sorted, unique and inside the interval
    //! \param[in] x The vector of points $x_i$ where to compute $s(x_i)$.
    //! \return Values of interpolant at $x$ (vector)
    Eigen::VectorXd operator() (Eigen::VectorXd x) const;

private:
    // Provided nodes and values $(t,y)$ to compute spline,
    // Eigen vectors, $c$ contains slopes
    // All have the same size
    Eigen::VectorXd t, y, c;
    // Difference $t(i)-t(i-1)$
    double h;
    // Size of $t$, $y$ and $c$.
    int n;
};
/* SAM_LISTING_BEGIN_0 */

PCHI::PCHI(const Eigen::VectorXd & t,
           const Eigen::VectorXd & y,
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
        /// CASE: assuming f'(x_j) = 0   (error: O(1))
        case Slope::Zero:
        c = Eigen::VectorXd::Zero(n);
            break;
    /* SAM_LISTING_END_1 */
    /* SAM_LISTING_BEGIN_2 */
        /// CASE: second order finite difference (error: O(h^2))
        case Slope::Reconstructed:
        default:

            // First order alternative:
//                c(0) = ( y(1) - y(0) ) / h;
            c(0) = ( -1*y(2) + 4*y(1) - 3*y(0) ) / 2 / h;
            for(int i = 1; i < n-1; ++i) {
                c(i) = ( y(i+1) - y(i-1) ) / 2 / h;
            }

            // First order alternative:
//                c(n-1) = ( y(n-1) - y(n-2) ) / h;
            c(n-1) = ( 3*y(n-1) - 4*y(n-2) + 1*y(n-3) ) / 2 / h;

            break;
    /* SAM_LISTING_END_2 */
    }

}

/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd PCHI::operator() (Eigen::VectorXd x) const {

    Eigen::VectorXd ret(x.size());

    for(int j = 0; j < x.size(); ++j) {
        // Stores the current interval index
        // and some temporary variable
        int i_star = 0;
        double tmp,t1,y1,y2,c1,c2,a1,a2,a3;
        // Find interval (porting of hermloceval Matlab code 3.4.6)
        if( t(i_star) < x(j) || i_star == 0) {
            ++i_star;
            t1 = t(i_star-1);
            y1 = y(i_star-1);
            y2 = y(i_star);
            c1 = c(i_star-1);
            c2 = c(i_star);
            a1 = y2 - y1;
            a2 = a1 - h*c1;
            a3 = h*c2 - a1 - a2;
        }
        // Compute $s(x(j))$
        tmp = ( x(j) - t1 ) / h;
        ret(j) = y1 + ( a1 + ( a2+a3*tmp ) * ( tmp - 1. ) ) * tmp;
    }

    return ret;
}
/* SAM_LISTING_BEGIN_2 */

int main() {
    // Interpoland
    auto f = [] (double x) { return 1. / (1. + x*x); };
    // auto f = [] (double x) {return cos(x); };

    double a = 5; // Interval bounds will be (-a,a)
    int M = 1000; // Number of  points in which to evaluate the interpoland

    // Precompute values at which evaluate f
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(M, -a, a);
    Eigen::VectorXd fx = x.unaryExpr(f);

    // Store error and number of nodes
    std::vector<double> N, err_reconstr, err_zero;

    for(int i = 2; i <= 512; i = i << 1) {
        // Define subintervals and evaluate f there (find pairs (t,y))
        Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(i, -a, a);
        Eigen::VectorXd y = t.unaryExpr(f);

        // Construct PCHI with zero and reconstructed slopes
        PCHI P_reconstr(t,y),
             P_zero(t, y, Slope::Zero);

        // Compute infinity norm of error
        Eigen::VectorXd P_zero_x = P_zero(x);
        Eigen::VectorXd P_reconstr_x = P_reconstr(x);
        err_reconstr.push_back(
                    (P_reconstr_x - fx).lpNorm<Eigen::Infinity>()
                    );
        err_zero.push_back(
                    (P_zero_x - fx).lpNorm<Eigen::Infinity>()
                    );
        N.push_back(i);

        // Se how interpolant looks
        if( i == 512 ) {
            {
                mgl::Figure fig;
                fig.title("Interpolant with zero slope");
                fig.xlabel("t");
                fig.ylabel("y");
                fig.plot(x, P_zero_x, "r").label("P_{zero}");
                fig.plot(x, fx, "b-").label("f");
                fig.legend();
                fig.save("p_zero");
            }
            {
                mgl::Figure fig;
                fig.title("Interpolant with reconstructed slope");
                fig.xlabel("t");
                fig.ylabel("y");
                fig.plot(x, P_reconstr_x, "r").label("P_{reconstr}");
                fig.plot(x, fx, "b-").label("f");
                fig.legend();
                fig.save("p_reconstr");
            }
        }
    }

    mgl::Figure fig;
    fig.title("Error VS no. of nodes");
    fig.setlog(true, true);
    fig.xlabel("No. of interpolation nodes");
    fig.ylabel("max |f(t) - P_\cdot(t)|");
    fig.plot(N, err_zero, "r").label("P_{zero}");
    fig.plot(N, err_reconstr, "b").label("P_{reconstr}");
    fig.legend();
    fig.save("pchi_conv");
}
