/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2022
 */

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace PolyEvalNewton {

// Class supporting Lagrange polynomial interpolation
/* SAM_LISTING_BEGIN_1 */
class PolyEval {
 public:
  // Idle Constructor
  PolyEval() = default;
  // Add another data point and update internal information
  void addPoint(double t, double y, double tol = 1.0E-10);
  // Evaluation of current interpolating polynomial at many points
  double operator()(double x) const;

 private:
  std::vector<double> t_;
  std::vector<double> c_;
};
/* SAM_LISTING_END_1 */

// Updating Newton representation of interpolating polynomial using the
// divided-difference algorithm
/* SAM_LISTING_BEGIN_2 */
void PolyEval::addPoint(double t, double y, double tol) {
  // Number of points already added
  const int n_pts = c_.size();
  // Store coordinates of new point
  t_.push_back(t);
  c_.push_back(y);
  // Update coefficients by means of divided-difference formula
  for (int j = n_pts - 1; j >= 0; --j) {
    const double td = t_[j] - t;
    // Make sure that the data abscissas are pairwise distinct
    if (std::abs(td) < tol * std::max(std::abs(t_[j]), std::abs(t))) {
      throw std::runtime_error("Duplicate data coordinate");
    }
    // Update \prbeqref{eq:ddfpr}
    c_[j] = (c_[j] - c_[j + 1]) / td;
  }
}
/* SAM_LISTING_END_2 */

// Evaluation operator implementing Horner scheme
/* SAM_LISTING_BEGIN_3 */
double PolyEval::operator()(double x) const {
  const int n = c_.size();
  double p = (n == 0) ? 0.0 : c_[0];
  for (int j = 1; j < n; ++j) {
    p = (x - t_[j]) * p + c_[j];
  }
  return p;
}
/* SAM_LISTING_END_3 */

// The same as a stand-alone function
/* SAM_LISTING_BEGIN_4 */
double evalNewtonForm(const std::vector<double> &t,
                      const std::vector<double> &c, double x) {
  assert((t.size() == c.size() && "Sizes of t and c vectors must agree"));
  const int n = c.size() - 1;
  double p = c[n];
  for (int j = n - 1; j >= 0; --j) {
    p = (x - t[j]) * p + c[j];
  }
  return p;
}
/* SAM_LISTING_END_4 */

// Generator class for sequence
class ConvergentSequence {
 public:
  ConvergentSequence() = default;
  double next() { return std::atan(n++); }
  void reset() { n = 1; }

 private:
  unsigned int n{1};
};

// Extrapolation to zero based on class PolyEval
/* SAM_LISTING_BEGIN_6 */
template <class SEQGEN>
double extrapolateLimit(SEQGEN &&seq, double rtol = 1.0E-6,
                        double atol = 1.0E-8, unsigned int maxdeg = 20) {
  seq.reset();
  PolyEval intpoly;  // Interpolating poklynomial
  intpoly.addPoint(1.0, seq.next());
  double lim_high = intpoly(0.0);
  for (unsigned i = 1; i < maxdeg; ++i) {
    const double lim_low = lim_high;
    // Add next data point $\cob{1/(i+1),s_{i+1}}$.
    intpoly.addPoint(1.0 / (1.0 + i), seq.next());
    // Evaluate polynomial interpolant in $\cob{x=0}$.
    lim_high = intpoly(0.0);
    // \com{error indicator} by comparing tow successive approximations
    const double errest = std::abs(lim_high - lim_low);
    if (errest < rtol * std::abs(lim_high) || errest < atol)  // \label{de:1}
      return lim_high;
  }
  // Specified tolerance could not be reached.
  throw std::runtime_error("extrapolateLimit failed to converge!");
  return std::nan("");
}
/* SAM_LISTING_END_6 */

// Extrapolation to zero for the computation of the limit of a sequence
/* SAM_LISTING_BEGIN_5 */
template <class SEQGEN>
double limitByExtrapolation(SEQGEN &&seq, double rtol = 1.0E-6,
                            double atol = 1.0E-8, unsigned int maxdeg = 20) {
  seq.reset();
  std::vector<double> h{1};
  std::vector<double> s{seq.next()};
  // using \com{Aitken-Neville scheme} with \Blue{$x=0$}, see
  // ~\lref{AitkenNeville}
  for (unsigned i = 1; i < maxdeg; ++i) {
    // Next data point $\cob{(h_i,s_i)}$
    h.push_back(1.0 / (i + 1.0));  // First coordinate of data point
    s.push_back(seq.next());       // Second coordinate of data point
    // Aitken-Neville update
    for (int k = i - 1; k >= 0; --k)
      s[k] = s[k + 1] - (s[k + 1] - s[k]) * h[i] / (h[i] - h[k]);
    // termination of extrapolation when desired tolerance is reached
    const double errest = std::abs(s[1] - s[0]);  // \com{error indicator}
    if (errest < rtol * std::abs(s[0]) || errest < atol) {  // \label{de:1}
      // Return value extrapolated from largest number of data points
      return s[0];
    }
  }
  throw std::runtime_error("limitByExtrapolation failed to converge!");
  return std::nan("");
}
/* SAM_LISTING_END_5 */

}  // namespace PolyEvalNewton

int main(int /*argc*/, char ** /*argv*/) {
  PolyEvalNewton::PolyEval polynom;
  polynom.addPoint(0.0, 1.0);
  polynom.addPoint(1.0, 2.0);
  polynom.addPoint(2.0, 5.0);
  std::cout << "p(1.5) = " << polynom(1.5) << std::endl;

  std::cout << "Limit = "
            << PolyEvalNewton::limitByExtrapolation(
                   PolyEvalNewton::ConvergentSequence())
            << std::endl;
  std::cout << "Limit = "
            << PolyEvalNewton::extrapolateLimit(
                   PolyEvalNewton::ConvergentSequence())
            << std::endl;
  return 0;
}
