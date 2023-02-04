/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2022
 * HW problem "On Quadrature Formulas" prb:qforms
 */

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace QuadFormulas {

/** Quadrature formula on reference interval [-1,1] */
/* SAM_LISTING_BEGIN_0 */
struct QuadFormulaRef {
  [[nodiscard]] virtual const std::vector<double>& nodes() const = 0;
  [[nodiscard]] virtual const std::vector<double>& weights() const = 0;
  virtual ~QuadFormulaRef() = default;
};
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
class QuadFormula {
 public:
  QuadFormula(const QuadFormulaRef& qf, double a, double b);
  QuadFormula(const QuadFormula&) = default;
  QuadFormula(QuadFormula&&) = default;
  QuadFormula& operator=(const QuadFormula&) = default;
  QuadFormula& operator=(QuadFormula&&) = default;

  [[nodiscard]] std::pair<double, double> interval() const { return {a_, b_}; }
  [[nodiscard]] const std::vector<double>& nodes() const { return nodes_; }
  [[nodiscard]] const std::vector<double>& weights() const { return weights_; }

 private:
  double a_;
  double b_;
  std::vector<double> nodes_;
  std::vector<double> weights_;
};
/* SAM_LISTING_END_1 */

// Constructor implementing the transformation of a quadrature formula
/* SAM_LISTING_BEGIN_2 */
QuadFormula::QuadFormula(const QuadFormulaRef& qf, double a, double b)
    : a_(a), b_(b), nodes_(qf.nodes()), weights_(qf.weights()) {
  const double len = b - a;
  const unsigned int n_pts = nodes_.size();
  assert((weights_.size() == n_pts) && "Size mismatch");
  for (unsigned int j = 0; j < n_pts; ++j) {
    // Affine transformation of quadrature nodes, see \lref{rem:quadtrf}.
    nodes_[j] = a + (nodes_[j] + 1.0) * len / 2;
    // Rescaling of weights
    weights_[j] *= len / 2;
  }
}
/* SAM_LISTING_END_2 */

// Determine maximal degree of polynomial integrated exactly
/* SAM_LISTING_BEGIN_3 */
unsigned int checkQuadOrder(const QuadFormulaRef& qf,
                            const double tol = 1.0E-10) {
  const std::vector<double>& c = qf.nodes();
  const std::vector<double>& w = qf.weights();
  const unsigned int n_pts = c.size();
  assert((w.size() == n_pts) && "Size mismatch");

  // Values of monomials at quadrature points
  std::vector<double> monom_vals(n_pts, 1.0);
  // Check for exactness for polynomials of degree d
  // Note that the highest possible order is $\cob{2N}$!
  for (unsigned int d = 0; d < n_pts; ++d) {
    // j-loop performs the evaluation of the quadrature formula
    // for the monomial $\cob{t \mapsto t^{2d}}$ (even polynomial)
    double s = 0.0;
    for (unsigned int j = 0; j < n_pts; ++j) {
      s += w[j] * monom_vals[j];
    }
    // Exact value: $\cob{2/(d+1)}$ for even degree, zero else.
    double val = 2.0 / (2 * d + 1);
    // Safe test for equality of computed floating point numbers
    if (std::abs(s - val) > tol) {
      return 2 * d;
    }
    // Now test for odd-degree monomials, for which the integral has to
    // evaluate to zero
    for (unsigned int j = 0; j < n_pts; ++j) {
      monom_vals[j] *= c[j];
    }
    s = 0.0;
    for (unsigned int j = 0; j < n_pts; ++j) {
      s += w[j] * monom_vals[j];
    }
    // Safe test for equality with zero
    if (std::abs(s) > tol) {
      return 2 * d + 1;
    }
    for (unsigned int j = 0; j < n_pts; ++j) {
      monom_vals[j] *= c[j];
    }
  }
  // The quadrature rule has the highest possible order
  return 2 * n_pts;
}
/* SAM_LISTING_END_3 */

class TwoPtGaussQF : public QuadFormulaRef {
 public:
  TwoPtGaussQF() {}
  TwoPtGaussQF(const TwoPtGaussQF&) = default;
  TwoPtGaussQF(TwoPtGaussQF&&) = default;
  TwoPtGaussQF& operator=(const TwoPtGaussQF&) = default;
  TwoPtGaussQF& operator=(TwoPtGaussQF&&) = default;
  virtual ~TwoPtGaussQF() = default;

  [[nodiscard]] virtual const std::vector<double>& nodes() const override {
    return nodes_;
  }
  [[nodiscard]] virtual const std::vector<double>& weights() const override {
    return weights_;
  }

 private:
  std::vector<double> nodes_{-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
  std::vector<double> weights_{1.0, 1.0};
};

}  // namespace QuadFormulas

int main(int /*argc*/, char** /*argv*/) {
  QuadFormulas::TwoPtGaussQF qf{};
  std::cout << "Order of 2-pt Gauss rule = " << QuadFormulas::checkQuadOrder(qf)
            << std::endl;

  return 0;
}
