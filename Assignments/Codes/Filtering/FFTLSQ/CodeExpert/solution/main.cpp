#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "fftlsq.hpp"
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

int main(int argc, char** argv) {
  // testing testNormEqMatrix
  unsigned int n = 10;
  unsigned int m = 3;
  bool test = testNormEqMatrix(n, m);
  if (test) {
    std::cout << "testNormEqMatrix passed!\n\n";
  } else {
    std::cout << "testNormEqMatrix failed!\n\n";
  }

  // Test points
  constexpr unsigned int npoints = 10;
  VectorXd d(npoints);
  d << 0.987214, 1.03579, 0.997689, 0.917471, 1.00474, 0.92209, 1.03517,
      1.08863, 0.904992, 0.956089;

  // Find best polynomial  (coefficients)
  VectorXd g;
  g = find_c(d, m);
  std::cout << "testing find_c" << std::endl;
  std::cout << g << std::endl << std::endl;

  plt::figure();
  // TODO: (5-2.e) Tabulate the coefficients of find_c for m = 1, 2, 3,
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
  plt::savefig("cx_out/orbit.eps");
}
