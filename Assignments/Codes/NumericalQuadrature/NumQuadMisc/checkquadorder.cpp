/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2021
 */

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/FFT>
#include <vector>

#define _USE_MATH_DEFINES

// Determining the order of a quadrature rule on [-1,1]
// by testing exact integration of monomials
/* SAM_LISTING_BEGIN_1 */
unsigned int checkQuadOrder(const Eigen::VectorXd &c, const Eigen::VectorXd &w,
                            const double tol = 1.0E-10) {
  const unsigned long N = c.size();  // Number of quadrature points
  assert(N == w.size());             // Numbers of nodes and weights must agree
  // Check for exactness for polynomials of degree d
  // Note that the highest possible order is $\cob{2N}$!
  for (unsigned int d = 0; d < 2 * N; ++d) {
    // j-loop performs the evaluation of the quadrature formula
    // for the monomial $\cob{t \mapsto t^d}$.
    double s = 0.0;
    for (int j = 0; j < N; ++j) {
      s += w[j] * std::pow(c[j], d);
    }
    // Exact value: $\cob{2/(d+1)}$ for even degree, zero else.
    double val = (d % 2 == 0) ? 2.0 / (d + 1) : 0.0;
    // Safe test for equality of computed floating point numbers
    if (std::abs(s - val) > tol) {
      return d;
    }
  }
  // The quadrature rule has the highest possible order
  return 2 * N;
}
/* SAM_LISTING_END_1 */

// Equidistant N-point Trapezoidal rule
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
double trapezoidalRule(FUNCTOR &&f, double a, double b, unsigned int N) {
  assert(a < b);
  assert(N >= 2);
  // Spacing of the N equidistant quadrature nodes
  const double h = (b - a) / (N - 1);
  // Variable containing location of current quadrature node
  double t = a + h;
  // Left endpoint has weight 0.5*h
  double s = 0.5 * f(a);
  // Summing contributions of N-2 interior quadrature nodes with weight h
  for (int j = 1; j < N - 1; ++j) {
    s += f(t);
    t += h;
  }
  // Right endpoint has weight 0.5*h
  s += 0.5 * f(t);
  // Scaling with h is performed in the end. Can also be done when accumulating
  // value of quadrature formula.
  return h * s;
}
/* SAM_LISTING_END_2 */

int main(int /*argc*/, char ** /*argv*/) {
  // 2-point Gauss quadrature rule of order 4
  std::cout << "Checking order of 2-point Gauss rule" << std::endl;
  Eigen::VectorXd w = Eigen::Vector2d(1.0, 1.0);
  Eigen::VectorXd c = Eigen::Vector2d(-std::sqrt(3) / 3.0, std::sqrt(3) / 3.0);
  std::cout << "Order of 2-pt Gauss rule = " << checkQuadOrder(c, w)
            << std::endl;

  std::cout << "Applying trapezoidal rule" << std::endl;
  std::cout << "int(-1,1,x -> x^2) = "
            << trapezoidalRule([](double t) { return t * t; }, -1.0, 1.0, 100)
            << std::endl;
  return 0;
}
