#include <iostream>
#include <vector>
#include <iomanip>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

//! \brief Function $\phi$ used in the Exponential Euler
//! single step method for an autonomous ODE.
MatrixXd phim(const MatrixXd & Z) {
    int n = Z.cols();
    assert( n == Z.rows() && "Matrix must be square.");
    MatrixXd C(2*n,2*n);
    C << Z, MatrixXd::Identity(n,n),
         MatrixXd::Zero(n,2*n);
    return C.exp().block(0, n, n, n);
}

//! \brief Calculate a single step of the exponential Euler method.
//! \tparam Function function object for r.h.s. function
//! \tparam Jacobian function object for Jacobian of r.h.s.
//! \param[in] y0 The initial state
//! \param[in] f The r.h.s function $f$
//! \param[in] df The Jacobian of $f$
//! \param[in] h The stepsize of the method.
//! \return A single step of the Exponential Euler method
/* SAM_LISTING_BEGIN_1 */
template <class Function, class Jacobian>
VectorXd ExpEulStep(const VectorXd & y0,
                    const Function& f, const Jacobian & df,
                    double h) {
   return  y0 + h * phim( h * df(y0)) * f(y0);
}
/* SAM_LISTING_END_1 */

int main() {
    /* SAM_LISTING_BEGIN_2 */
    // Final time
    double T = 1;
    // Initial value
    VectorXd y0(1);
    y0 << 0.1;
    // Function and Jacobian and exact solution
    auto f = [] (const VectorXd & y) {
        return y(0) * ( 1 - y(0) );
    };
    auto df = [] (const VectorXd & y) {
        VectorXd dfy(1);
        dfy << 1 - 2*y(0);
        return dfy;
    };
    double exactyT = y0(0) / (y0(0)+(1-y0(0))*exp(-T));

    // Container for errors
    std::vector<double> error(15);

    // Test many step sizes
    for (int j=0; j < 15; j++) {
        int N = std::pow(2,j+1);
        double h = T / N;
        VectorXd y = y0;
        for (int k=0; k < N; k++) {
            y = ExpEulStep(y,f,df,h);
        }

        error[j] = std::abs(y(0) - exactyT);
        std::cout << std::left << std::setfill(' ')
                  << std::setw(3) << "N = "
                  << std::setw(7) << N
                  << std::setw(8) << "Error = "
                  << std::setw(13) << error[j];
        if (j > 0) {
            std::cout << std::left << std::setfill(' ')
                      << std::setw(10) << "Approximated order = "
                      << std::log2( error[j-1]/error[j] )
                      << std::endl;
        }
        else std::cout << std::endl;
    }
    /* SAM_LISTING_END_2 */
}
