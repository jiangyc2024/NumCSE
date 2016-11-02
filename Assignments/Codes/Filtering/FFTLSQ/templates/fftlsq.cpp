#include <iostream>

#include <Eigen/Dense>

#include "figure/figure.hpp"
#include "FFT/fft.hpp"

using namespace Eigen;

/*!
 * \brief gauss_fit Find best trigonometric polynomial
 * passing trough distances $\mathbf{d}$.
 * Using Least Squares formulation and assuming best fitting
 * polynomial has degree $m$.
 * \param d Vector of size $n$ distances at angle $2*\pi*i/n$.
 * \param m Degree of the trigonometric polynomial $p$.
 * \return The coefficients of the trigonometric polinomial.
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd gauss_fit(const VectorXd & d,
                   unsigned int m) {
    unsigned int n = d.size();

    // We will use a real to complex, discrete Fourier transform.
    FFT<double> fft;
    VectorXd rhs;
    // TODO: implement this function
    return rhs;

    return rhs;
}
/* SAM_LISTING_END_1 */


/*!
 * \brief eval_p Given polynomial coefficients, return value of polynomial
 * at $n$ equidistant points.
 *
 * \param p Coefficient vector of trigonometrix polynomial.
 * \param n Number of equidistant points at which to evaluate.
 * \return Value of polynomial $p$ at $2\pi i / n$.
 */
VectorXd eval_p(VectorXd c, unsigned int n) {

    // Degree of polynomial
    unsigned int m = c.size();

    VectorXd ret(n);
    // Loop over all points
    for (unsigned int i = 0; i < n; ++i) {
        double r = 0;
        // Loop over all coefficients
        for (unsigned int j = 0; j < m; ++j) {
            r += c(j) * std::cos(2 * M_PI * i * j / n);
        }
        ret(i) += r;
    }
    return ret;
}

int main(int argc, char **argv) {

    // Degree of trigonometric polynomial
    unsigned int m = 3;
    if(argc > 1) {
        m = std::stoi(argv[1]);
    }

    // Test points
    unsigned int npoints = 10;
    VectorXd d(npoints);
    d << 0.987214,
         1.03579,
        0.997689,
        0.917471,
         1.00474,
         0.92209,
         1.03517,
         1.08863,
        0.904992,
        0.956089;

    VectorXd g;
    // TODO: Find coefficients that best fit the data in d
    // using trig. poly of degree $m$, store the result in $g$.

    // Find coordinates of best poly coeff.
    unsigned int neval = 100;
    VectorXd e = eval_p(g, neval);
    VectorXd x, y;
    x.resizeLike(e);
    y.resizeLike(e);
    for(unsigned int i = 0; i < neval; ++i) {
        x(i) = std::sin(2. * M_PI * i / neval) * e(i);
        y(i) = std::cos(2. * M_PI * i / neval) * e(i);
    }

    // Find coordinates of points
    VectorXd x_p, y_p;
    x_p.resizeLike(d);
    y_p.resizeLike(d);
    for(unsigned int i = 0; i < npoints; ++i) {
        x_p(i) = std::sin(2 * M_PI * i / npoints) * d(i);
        y_p(i) = std::cos(2 * M_PI * i / npoints) * d(i);
    }

    // Plot points and poly
    mgl::Figure fig;
    fig.title("Orbit of planet");
    fig.ranges(-2, 2, -2, 2);
    fig.plot(x, y, "r").label("best orbit");
    fig.plot(x_p, y_p, " b*").label("points");
    fig.xlabel("x");
    fig.ylabel("y");
    fig.legend(0, 1);
    fig.save("orbit.eps");
    fig.save("orbit.png");

}
