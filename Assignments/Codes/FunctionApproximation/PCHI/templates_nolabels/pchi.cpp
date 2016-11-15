//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <vector>
#include <cassert>

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <figure.hpp>

// Flag for slope reconstruction type
enum class Slope { Zero, Reconstructed };

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
        /// CASE: assuming $s'(x_j) = 0$ (error: $O(1)$)
        case Slope::Zero:
    // TODO: set the slopes assuming $s'(x_j) = 0$
            break;
        /// CASE: second order finite differences (error: $O(h^2)$)
        case Slope::Reconstructed:
        default:

        // TODO: set the slopes using finite difference reconstruction
            break;
    }

}

Eigen::VectorXd PCHI::operator() (Eigen::VectorXd x) const {

    Eigen::VectorXd ret(x.size());



    return ret;
}

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
