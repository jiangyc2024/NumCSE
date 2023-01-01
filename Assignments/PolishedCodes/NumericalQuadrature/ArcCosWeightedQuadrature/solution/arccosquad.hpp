#ifndef ARCCOSQUAD_HPP
#define ARCCOSQUAD_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// The function gaussquad() for pre-computed (and pre-compiled)
// Gauss-Legendre quadrature rules up to 256 nodes.
#include "gaussquad.hpp"

/* SAM_LISTING_BEGIN_0 */
struct cvgdata {
  long n_{0};          // number of quadrature points
  double val_{0.0};    // approximate value
  double err_{-1.0};   // quadrature error
  double rate_{-1.0};  // estimate for rate of algebraic convergence
};
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void testConvGaussQuad() {
  // Table header
  std::cout << std::setw(5) << "n" << std::setw(20) << " I_n " << std::setw(30)
            << " Error " << std::setw(10) << (" rate\n " + std::string(65, '-'))
            << std::endl;

  // Value obtained by overkill quadrature as reference value
  constexpr double ref_val = 1.7576137811123187;

  // TODO: (7-13.a) Print a table that allows you to predict the asymptotic
  // behaviour of Gauss-Legendre numerical quadrature when approximating I(f).
  // START
  // The integrand as lambda funtion
  auto integrand = [](double t) {
    return std::acos(t) * 1.0 / (1.0 + std::exp(t));
  };

  // For recording errors etc.
  std::vector<cvgdata> data;
  // We suspect algebraic convergence. Accordingly use doubling of number of
  // quadrature nodes
  constexpr unsigned int L = 8;  // Use at most $2^L$ quadrature nodes
  unsigned int n = 2;
  for (unsigned int l = 1; l <= L; ++l, n *= 2) {
    const QuadRule qr{gaussquad(n)};
    // No transformation required: the integral is on [-1,1] anyway.
    double I = 0.0;  // Summation variable for quadrature formula
    for (unsigned int j = 0; j < n; ++j) {
      I += qr.weights_[j] * integrand(qr.nodes_[j]);
    }
    data.push_back({qr.weights_.size(), I, 0.0, 0.0});
  }
  for (unsigned int l = 0; l < data.size(); ++l) {
    data[l].err_ = std::abs(data[l].val_ - ref_val);
    if (l > 0) {
      // Estimate for rate of algebraic convergence
      data[l].rate_ =
          (std::log(data[l - 1].err_) - std::log(data[l].err_)) / std::log(2.0);
    }
  }

  // Output table
  for (unsigned int l = 0; l < data.size(); ++l) {
    const cvgdata& item{data[l]};
    std::cout << std::setw(5) << item.n_ << std::setw(20) << std::fixed
              << std::setprecision(16) << item.val_ << std::setw(30)
              << std::scientific << std::setprecision(16) << item.err_
              << std::setw(10) << std::fixed << std::setprecision(2)
              << item.rate_ << std::endl;
  }
  // END
}
/* SAM_LISTING_END_1 */

/**
 * \brief Calculates the quadrature of f using transformation.
 *
 * \tparam FUNCTION suitable functor
 * \param f function to be integrated on [-1, 1]
 * \param n number of evaluations
 * \return double integral
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTION>
double arccosWeightedQuad(FUNCTION&& f, unsigned int n) {
  double I = 0.0;  // For accumulating quadrature result
  // TODO: (7-13.c) Approximate I(f) with exponential convergence in n.
  // START
  // Obtain Gauss-Legendre quadrature rule
  QuadRule qr{gaussquad(n)};
  // The transformed integral covers the interval $\cob{[0,\pi]}$, which entails
  // a transformation of the quadrature nodes and weights, see
  // \lref{rem:quadtrf}
  qr.weights_ *= 0.5 * M_PI;
  qr.nodes_ = 0.5 * M_PI * (1.0 + qr.nodes_.array());
  // Straightforward implementation of quadrature formula
  for (unsigned int j = 0; j < n; ++j) {
    const double c = qr.nodes_[j];
    const double fval = f(std::cos(c));
    I += qr.weights_[j] * fval * std::sin(c) * c;
  }
  // END
  return I;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void testConvTrfGaussQuad() {
  // Table header
  std::cout << std::setw(5) << "n" << std::setw(20) << " I_n " << std::setw(30)
            << " Error " << std::setw(10)
            << (" decay ratio\n " + std::string(65, '-')) << std::endl;

  // Value obtained by overkill quadrature as reference value
  constexpr double ref_val = 1.7576137811123187;

  // TODO: (7-13.d) Print a table that allows you to predict the asymptotic
  // behaviour of arccosWeightedQuad().
  // START
  // The integrand (/arccos(t)) as lambda funtion
  auto f = [](double t) { return 1.0 / (1.0 + std::exp(t)); };

  // For recording errors etc.
  std::vector<cvgdata> data;
  // We expect exponential convergence. Increment number of quadrature nodes by
  // one in each step
  constexpr unsigned int L = 12;  // Use a most L quadrature nodes
  for (unsigned int n = 2; n <= L; ++n) {
    data.push_back({n, arccosWeightedQuad(f, n), 0.0, 0.0});
  }
  for (unsigned int l = 0; l < data.size(); ++l) {
    data[l].err_ = std::abs(data[l].val_ - ref_val);
    if (l > 0) {
      // Quotient of successive error values
      data[l].rate_ = data[l].err_ / data[l - 1].err_;
    }
  }

  // Output table
  for (unsigned int l = 0; l < data.size(); ++l) {
    const cvgdata& item{data[l]};
    std::cout << std::setw(5) << item.n_ << std::setw(20) << std::fixed
              << std::setprecision(16) << item.val_ << std::setw(30)
              << std::scientific << std::setprecision(16) << item.err_
              << std::setw(10) << std::fixed << std::setprecision(2)
              << item.rate_ << std::endl;
  }
  // END
}
/* SAM_LISTING_END_3 */

#endif