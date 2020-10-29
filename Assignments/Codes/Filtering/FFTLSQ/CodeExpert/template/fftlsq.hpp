#ifndef FFTLSQ_HPP
#define FFTLSQ_HPP

#include <math.h>

#include <Eigen/Dense>
#include <iomanip>
#include <unsupported/Eigen/FFT>
#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;

/*!
 * \brief eval_p Given polynomial coefficients, return value of polynomial
 * at $n$ equidistant points.
 *
 * \param p Coefficient vector of trigonometrix polynomial.
 * \param n Number of equidistant points at which to evaluate.
 * \return Value of polynomial $p$ at $2\pi i / n$.
 */
VectorXd eval_p(const VectorXd& c, const unsigned int n) {
  // Degree of polynomial
  const unsigned int m = c.size();

  VectorXd ret(n);
  // Loop over all points
  for (unsigned int i = 0; i < n; ++i) {
    double r = 0;
    // Loop over all coefficients
    for (unsigned int j = 0; j < m; ++j) {
      r += c(j) * std::cos(2 * M_PI * i * j / n);
    }
    ret(i) = r;
  }
  return ret;
}


/*!
 * \brief testNormEqMatrix Create the matrix $A^TA$
 * in two different ways and make sure they are
 * approximately equal.
 * \param n number of different measurements
 * \param m degree of the trigonometric polynomial to
 * be fitted
 */
/* SAM_LISTING_BEGIN_0 */
bool testNormEqMatrix(unsigned int n, unsigned int m) {
  // TODO: (5-2.c) Test if the two definitions, given in (5.2.4) and
  // (5.2.7), of $A^T*A$ are equal.
  // START
  MatrixXd A(n, m + 1);
  // Initializing the vectors described in (5.2.4)
  VectorXd n_linear = VectorXd::LinSpaced(n, 0, n - 1);
  VectorXd m_linear = VectorXd::LinSpaced(m + 1, 0, m);
  // building A as described in (5.2.4)
  A = 2.0 * M_PI / n * n_linear * m_linear.transpose();
  A = A.array().cos().matrix();
  // Using the knowledge gained in (5-2.b) to
  // initialize $A^T*A$ directly as a diagonal matrix.
  MatrixXd ATA = MatrixXd::Zero(m + 1, m + 1);
  ATA.diagonal() = VectorXd::Constant(m + 1, n / 2.0);
  ATA(0, 0) = n;
  // Comparing the two approaches to calculate $A^T*A$.
  if ((ATA - A.transpose() * A).norm() < 1e-10 * m * n) return true;
  // END
  return false;
}
/* SAM_LISTING_0 */

/*!
 * \brief find_c Find best trigonometric polynomial
 * passing trough distances $\mathbf{d}$.
 * Using Least Squares formulation and assuming best fitting
 * polynomial has degree $m$.
 * \param d Vector of size $n$ distances at angle $2*\pi*i/n$.
 * \param m Degree of the trigonometric polynomial $p$.
 * \return The coefficients of the trigonometric polynomial.
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd find_c(const VectorXd& d, unsigned int m) {
  unsigned int n = d.size();

  // We will use a real to complex, discrete Fourier transform.
  FFT<double> fft;
  VectorXd rhs;
  // TODO: (5-2.d) find the coefficients
  // START
  // Computing DFT of d
  VectorXcd fourier = fft.fwd(d);
  // Gathering what is important for the rhs.
  rhs = fourier.real().head(m);
  // Solving normal equation by inverting diagonal matrix
  rhs /= n;
  rhs.tail(m - 1) *= 2;
  // END
  return rhs;
}
/* SAM_LISTING_END_1 */


/*!
 * \brief using implementation of find_c compute coefficients c_k and p* for m=1,2,3 
 * using matplotlibcpp's plt() create a plot showing the ellipse
 * in the same plot insert the curves descriped by p* of degrees m=1,2,3
 */
/* SAM_LISTING_BEGIN_2 */
void fitEllipse(void) {
    
  unsigned int n = 10;
  unsigned int m = 3;
  
  // Test points
  constexpr unsigned int npoints = 10;
  VectorXd d(npoints);
  d << 0.987214, 1.03579, 0.997689, 0.917471, 1.00474, 0.92209, 1.03517,
      1.08863, 0.904992, 0.956089;
  

  plt::figure();
  // TODO: (5-2.e) Tabulate the coefficients of find\_c for all $m=1,2,3$
  // plot the ellipse and also the curves of the trigonimetric polynomials
  // START
  constexpr double c = 0.8;
  auto d_func =
      [](const double phi) {  // no need to capture c because it is constexpr
        assert(0. <= phi && phi <= 2. * M_PI);
        return 1. / std::sqrt(1. - std::pow(c * std::cos(phi), 2));
      };

  // sample distance data, we reuse d from above
  n = 100;
  d.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    d(i) = d_func((2. * M_PI * i) / n);
  }

  plt::title("Orbit of planet");
  plt::xlim(-2, 2);
  plt::ylim(-2, 2);

  // plot the ellipse given by (5.2.10) using 10000 points
  constexpr unsigned int N = 10000;
  VectorXd x_ellipse(N), y_ellipse(N);
  double temp_d;
  for (std::size_t i = 0; i < N; ++i) {
    temp_d = d_func((2. * M_PI * i) / N);
    x_ellipse(i) = temp_d * std::cos((2 * M_PI * i) / N);
    y_ellipse(i) = temp_d * std::sin((2 * M_PI * i) / N);
  }
  plt::plot(x_ellipse, y_ellipse, "k", {{"label", "exact orbit"}});

  // calculate, tabulate and plot the coefficients in one go
  x_ellipse.resize(n);
  y_ellipse.resize(n);
  // create vector of colors for plots
  std::vector<std::string> colors = {"m", "y--", "r"};
  std::cout << "degree of polynomial | coefficients" << std::endl;
  for (m = 1; m <= 3; ++m) {
    VectorXd coefficients = find_c(d, m);

    std::cout << std::setw(20) << m << " | " << std::setw(5)
              << coefficients.transpose() << std::endl;

    VectorXd e = eval_p(coefficients, n);
    for (std::size_t i = 0; i < n; ++i) {
      x_ellipse(i) = e(i) * std::cos((2 * M_PI * i) / n);
      y_ellipse(i) = e(i) * std::sin((2 * M_PI * i) / n);
    }
    const std::string label = "fitted trig. poly. for m = " + std::to_string(m);
    plt::plot(x_ellipse, y_ellipse, colors[m - 1], {{"label", label}});
  }

  plt::xlabel("x");
  plt::ylabel("y");
  plt::legend();
  // END
  plt::savefig("cx_out/orbit.png");

}

/* SAM_LISTING_END_2 */

#endif
