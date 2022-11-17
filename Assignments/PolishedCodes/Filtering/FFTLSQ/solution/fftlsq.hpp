#ifndef FFTLSQ_HPP
#define FFTLSQ_HPP

#include <math.h>

#include <Eigen/Dense>
#include <iomanip>
#include <unsupported/Eigen/FFT>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * @brief eval_p Given polynomial coefficients, return value of polynomial
 * at $n$ equidistant points.
 *
 * @param c Coefficient vector of trigonometrix polynomial.
 * @param n Number of equidistant points at which to evaluate.
 * @return Eigen::VectorXd Value of polynomial $p$ at $2\pi i / n$.
 */
Eigen::VectorXd eval_p(const Eigen::VectorXd &c, const unsigned int n) {
  // Degree of polynomial
  const unsigned int m = c.size();

  Eigen::VectorXd ret(n);
  // Loop over all points
  for (unsigned int j = 0; j < n; ++j) {
    double r = 0;
    // Loop over all coefficients
    for (unsigned int k = 0; k < m; ++k) {
      r += c(k) * std::cos(2.0 * M_PI * k * j / n);
    }
    ret(j) = r;
  }

  return ret;
}

/**
 * @brief testNormEqMatrix Create the matrix $A^TA$ in two different ways and
 * make sure they are approximately equal.
 *
 * @param n number of different measurements
 * @param m degree of the trigonometric polynomial to
 * be fitted
 * @return bool indicating if test passed or failed.
 */
/* SAM_LISTING_BEGIN_0 */
bool testNormEqMatrix(unsigned int n, unsigned int m) {
  // TODO: (4-2.c) Test if the two definitions, given in (4.2.4) and
  // (4.2.7), of $A^T*A$ are equal.

  // START
  Eigen::MatrixXd A(n, m + 1);

  // Initializing the vectors described in (4.2.4)
  Eigen::VectorXd n_linear = Eigen::VectorXd::LinSpaced(n, 0, n - 1);
  Eigen::VectorXd m_linear = Eigen::VectorXd::LinSpaced(m + 1, 0, m);

  // Building A as described in (4.2.4)
  A = 2.0 * M_PI / n * n_linear * m_linear.transpose();
  A = A.array().cos().matrix();

  // Using the knowledge gained in (4-2.b) to initialize $A^T*A$ directly as a
  // diagonal matrix.
  Eigen::MatrixXd ATA = Eigen::MatrixXd::Zero(m + 1, m + 1);
  ATA.diagonal() = Eigen::VectorXd::Constant(m + 1, n / 2.0);
  ATA(0, 0) = n;

  // Comparing the two approaches to calculate $A^T*A$.
  // Test for relative smallness of the Frobenius norm of the difference of the
  // matrices
  if ((ATA - A.transpose() * A).norm() < 1e-10 * m * n) {
    return true;
  }
  // END

  return false;
}
/* SAM_LISTING_END_0 */

/**
 * @brief find_c Find Best trigonometric polynomial passing trough distances
 * $\mathbf{d}$. Using Least Squares formulation and assuming best fitting
 * polynomial has degree $m$.
 *
 * @param d Vector of size $n$ distances at angle $2*\pi*i/n$.
 * @param m Degree of the trigonometric polynomial $p$.
 * @return Eigen::VectorXd The coefficients of the trigonometric polynomial.
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd find_c(const Eigen::VectorXd &d, unsigned int m) {
  unsigned int n = d.size();

  // We will use a real to complex, discrete Fourier transform.
  Eigen::FFT<double> fft;
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m + 1);

  // TODO: (4-2.d) find the coefficients
  // START

  // Computing DFT of d
  Eigen::VectorXcd fourier = fft.fwd(d);

  // Gathering what is important for the rhs.
  rhs = fourier.real().head(m + 1);

  // Solving normal equation by inverting diagonal matrix.
  rhs /= n;
  rhs.tail(m) *= 2.0;
  // END

  return rhs;
}
/* SAM_LISTING_END_1 */

/**
 * @brief using implementation of find_c compute coefficients c_k and p* for
 * m=1,2,3 using matplotlibcpp's plt() create a plot showing the ellipse in the
 * same plot insert the curves descriped by p* of degrees m=1,2,3
 *
 */
/* SAM_LISTING_BEGIN_2 */
void fitEllipse() {
  unsigned int n = 10;
  unsigned int m = 3;

  // Test points
  constexpr unsigned int npoints = 10;
  Eigen::VectorXd d(npoints);
  d << 0.987214, 1.03579, 0.997689, 0.917471, 1.00474, 0.92209, 1.03517,
      1.08863, 0.904992, 0.956089;

  plt::figure();

  // TODO: (4-2.e) Tabulate the coefficients of find\_c for all $m=1,2,3$
  // plot the ellipse and also the curves of the trigonimetric polynomials
  // START
  constexpr double c = 0.8;
  auto d_func =
      [](const double phi) {  // no need to capture c because it is constexpr
        assert(0. <= phi && phi <= 2. * M_PI);
        return 1. / std::sqrt(1. - std::pow(c * std::cos(phi), 2));
      };

  // Sample distance data, we reuse d from above
  n = 100;
  d.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    d(i) = d_func((2. * M_PI * i) / n);
  }

  plt::title("Orbit of planet");
  plt::xlim(-2, 2);
  plt::ylim(-2, 2);

  // Plot the ellipse given by (5.2.10) using 10000 points
  constexpr unsigned int N = 10000;
  Eigen::VectorXd x_ellipse(N), y_ellipse(N);
  double temp_d;
  for (std::size_t i = 0; i < N; ++i) {
    temp_d = d_func((2.0 * M_PI * i) / N);
    x_ellipse(i) = temp_d * std::cos((2.0 * M_PI * i) / N);
    y_ellipse(i) = temp_d * std::sin((2.0 * M_PI * i) / N);
  }
  plt::plot(x_ellipse, y_ellipse, "k", {{"label", "exact orbit"}});

  // Calculate, tabulate and plot the coefficients in one go
  x_ellipse.resize(n);
  y_ellipse.resize(n);

  // Create vector of colors for plots
  std::vector<std::string> colors = {"m", "y--", "r"};
  std::cout << "degree of polynomial | coefficients" << std::endl;
  for (m = 1; m <= 3; ++m) {
    Eigen::VectorXd coefficients = find_c(d, m);

    std::cout << std::setw(20) << m << " | " << std::setw(5)
              << coefficients.transpose() << std::endl;

    Eigen::VectorXd e = eval_p(coefficients, n);
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
