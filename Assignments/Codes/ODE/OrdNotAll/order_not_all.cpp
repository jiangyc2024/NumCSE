#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>

#include "rkintegrator.hpp"

using namespace Eigen;

/*!
 * \brief errors Approximates the order of convergence of a scheme.
 * The scheme is a Runge Kutta scheme defined by A and b when applied
 * to the first order system y'=f(y), y(0)=y0.
 * The function ouputs error of the solutions at the time T.
 * \tparam Function Type for r.h.s function f.
 * \param f
 * \param T
 * \param y0
 * \param A
 * \param b
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function>
void errors(const Function &f, double T,
            const VectorXd &y0,
            const MatrixXd &A, const VectorXd &b) {

    RKIntegrator<VectorXd> rk(A, b);

    std::vector<double> error(15);
#if SOLUTION
    std::vector<double> order(14);

    double sum = 0;
    int count = 0;
    bool test = 1;

    std::vector<VectorXd> y_exact = rk.solve(f, T, y0, pow(2,15));

    for(int k = 0; k < 15; k++) {
        int N = std::pow(2,k+1);
        std::vector<VectorXd> y1 = rk.solve(f,T,y0,N);

        error[k] = (y1[N] - y_exact[std::pow(2,15)]).norm();

        std::cout << std::left << std::setfill(' ')
                  << std::setw(3)  << "N = "
                  << std::setw(7) << N
                  << std::setw(8) << "Error = "
                  << std::setw(13) << error[k];

        if ( error[k] < y0.size() * 5e-14) {
            test = 0;
        }
        if ( k>0 && test ) {
            order[k-1]=std::log2( error[k-1] / error[k] );
            std::cout << std::left << std::setfill(' ')
                      << std::setw(10)
                      << "Approximated order = " << order[k-1]
                      << std::endl;
            sum += order[k-1];
            count = k;
        }
        else std::cout << std::endl;
    }
    std::cout << "Average approximated order = "
              << sum / count
              << std::endl << std::endl;
#else // TEMPLATE
    // TODO: output error and order of the method
#endif
}
/* SAM_LISTING_END_1 */

int main() {
    /* SAM_LISTING_BEGIN_2 */
#if SOLUTION
    // Construct data for Butcher schemes
    MatrixXd A1 = MatrixXd::Zero(1,1);
    VectorXd b1(1);
    b1 << 1;

    MatrixXd A2 = MatrixXd::Zero(2,2);
    A2(1,0) = 1;
    VectorXd b2(2);
    b2 << .5, .5;

    MatrixXd A3 = MatrixXd::Zero(3,3);
    A3(1,0) = .5;
    A3(2,0) = -1;
    A3(2,1) = 2;
    VectorXd b3(3);
    b3 << 1./6, 2./3, 1./6;

    MatrixXd A4 = MatrixXd::Zero(4,4);
    A4(1,0) = .5;
    A4(2,1) = .5;
    A4(3,2) = 1;
    VectorXd b4(4);
    b4 << 1./6, 1./3, 1./3, 1./6;

    // First ODE
    cout << endl << "1. ODE y' = (1-y)y, y(0)=.5" << endl << endl;
    double T = 0.1;
    auto f = [] (VectorXd y) {VectorXd fy(1); fy << (1.-y(0))*y(0); return fy;};
    VectorXd y0(1); y0 << .5;

    cout << "Explicit Euler"
         << endl << endl;
    errors(f, T, y0, A1, b1);
    cout << "Trapezoidal rule"
         << endl << endl;
    errors(f, T, y0, A2, b2);
    cout << "RK order 3"
         << endl << endl;
    errors(f, T, y0, A3, b3);
    cout << "Classical RK order 4"
         << endl << endl;
    errors(f, T, y0, A4, b4);

    // Second ODE
    cout << endl << "2. ODE y' = |1.1 - y| + 1, y(0)=1" << endl << endl;
    auto f2 = [] (VectorXd y) {VectorXd fy(1); fy << abs(1.1-y(0))+1.; return fy;};
    y0 << 1;

    cout << "Explicit Euler"
         << endl << endl;
    errors(f2, T, y0, A1, b1);
    cout << "Trapezoidal rule"
         << endl << endl;
    errors(f2, T, y0, A2, b2);
    cout << "RK order 3"
         << endl << endl;
    errors(f2, T,  y0, A3, b3);
    cout << "Classical RK order 4"
         << endl << endl;
    errors(f2, T, y0, A4, b4);
#else // TEMPLATE
    // TODO: output error and order of the method
#endif
    /* SAM_LISTING_END_2 */
}
