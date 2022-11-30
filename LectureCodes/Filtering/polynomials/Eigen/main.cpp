/* Demonstration code for Lecture Numerical Methods for CSE
 * Polynomials arithmetic
 * Author: R. Hiptmair
 * Date: October 2020
 * (C) SAM, ETZ Zurich
 */

#include <Eigen/Dense>
#include <cassert>
#include <complex>
#include <iostream>
#include <unsupported/Eigen/FFT>
#include <utility>
#include <vector>

namespace Polynomials {
/* Note: Polynomials are encoded as vectors of coefficients in monomial
 * representation, the coefficient of the constant term providing the v[0] */

// Multiplication of two polynomials
template <typename RetVec, typename uVec, typename vVec>
RetVec polyMult(const uVec &u, const vVec &v) {
  // Fetch degrees of input polynomials
  const int degu = u.size() - 1;
  const int degv = v.size() - 1;
  // Object for product polynomial p = uv
  const int degp = degu + degv;
  RetVec p(degp + 1);
  // Perform convolution
  for (int j = 0; j <= degp; ++j) {
    const int l_low = std::max(0, j - degv);
    const int l_high = std::min(degu, j);
    for (int l = l_low; l <= l_high; ++l) {
      p[j] += u[l] * v[j - l];
    }
  }
  return p;
}

// Euclids algorithm: division of polynomials: p = u*x+r
template <typename RetVec, typename pVec, typename uVec>
std::pair<RetVec, RetVec> polyDiv(pVec p, const uVec &u) {
  using uScalar = typename uVec::value_type;
  using pScalar = typename pVec::value_type;
  // Fetch degrees of input polynomials
  const int degp = p.size() - 1;
  int degu = u.size() - 1;
  assert(degp >= degu && "Degree of p must be >= degree of u");
  // Effective degree of u
  while (u[degu] == static_cast<uScalar>(0)) {
    degu--;
  }
  assert(degu > 0 && "u != 0 required");
  // Degree of result of division
  const int degv = degp - degu;
  RetVec v(degv + 1);
  // Degree of remainder
  const int degr = degp - 1;
  RetVec r(degr + 1);
  // Euclid's algorithm
  for (int i = 0; i <= degv; ++i) {
    const pScalar sigma = p[degp - i] / u[degu];
    for (int l = 0; l <= degu; ++l) {
      p[degp - i - degu + l] -= sigma * u[l];
    }
    v[degv - i] = sigma;
  }
  for (int i = 0; i <= degr; ++i) {
    r[i] = p[i];
  }
  return {v, r};
}

// Print coefficients of a polynomial
template <typename Vec> void polyOut(const Vec &v) {
  std::cout << '[';
  for (auto x : v) {
    std::cout << x << ' ';
  }
  std::cout << ']' << std::endl;
}
} // namespace Polynomials

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "Polynomial multiplication and division" << std::endl;
  const std::vector<double> u{1.0, 2.0, 3.0};
  const std::vector<double> v{4.0, 5.0};
  auto p = Polynomials::polyMult<std::vector<double>>(u, v);
  std::cout << "p = ";
  Polynomials::polyOut(p);
  auto [q, r] = Polynomials::polyDiv<std::vector<double>>(p, u);
  std::cout << "q = ";
  Polynomials::polyOut(q);
  std::cout << "r = ";
  Polynomials::polyOut(r);
  return 0;
}
