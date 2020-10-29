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

  // END
  plt::savefig("cx_out/orbit.eps");
}
