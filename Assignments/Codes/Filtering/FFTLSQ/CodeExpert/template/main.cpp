#include <Eigen/Dense>
#include <iostream>

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
  // testing testNormEqMatrix
  unsigned int n = 10;
  unsigned int m = 3;
  bool test = testNormEqMatrix(n, m);
  if (test) {
    std::cout << "testNormEqMatrix passed!\n\n";
  } else {
    std::cout << "testNormEqMatrix failed!\n\n";
  }
  // Degree of trigonometric polynomial
  if (argc > 1) {
    m = std::stoi(argv[1]);
  }

  // Test points
  unsigned int npoints = 10;
  VectorXd d(npoints);
  d << 0.987214, 1.03579, 0.997689, 0.917471, 1.00474, 0.92209, 1.03517,
      1.08863, 0.904992, 0.956089;

  // Find best polynomial  (coefficients)
  VectorXd g;
  g = find_c(d, m);
  std::cout << "testing find_c" << std::endl;
  std::cout << g << std::endl << std::endl;

  // Find coordinates of best poly coeff.
  unsigned int neval = 100;
  VectorXd e = eval_p(g, neval);
  VectorXd x, y;
  x.resizeLike(e);
  y.resizeLike(e);
  for (unsigned int i = 0; i < neval; ++i) {
    x(i) = std::sin(2. * M_PI * i / neval) * e(i);
    y(i) = std::cos(2. * M_PI * i / neval) * e(i);
  }

  // Find coordinates of points
  VectorXd x_p, y_p;
  x_p.resizeLike(d);
  y_p.resizeLike(d);
  for (unsigned int i = 0; i < npoints; ++i) {
    x_p(i) = std::sin(2 * M_PI * i / npoints) * d(i);
    y_p(i) = std::cos(2 * M_PI * i / npoints) * d(i);
  }

  // Plot points and poly
  plt::figure();
  plt::title("Orbit of planet");
  plt::xlim(-2, 2);
  plt::ylim(-2, 2);
  plt::plot(x, y, "r", {{"label", "best orbit"}});
  plt::plot(x_p, y_p, " b*", {{"label", "points"}});
  plt::xlabel("x");
  plt::ylabel("y");
  plt::legend();
  plt::savefig("cx_out/orbit.png");
}
