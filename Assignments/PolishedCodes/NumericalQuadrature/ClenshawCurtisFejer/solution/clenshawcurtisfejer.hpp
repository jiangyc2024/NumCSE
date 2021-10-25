#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <unsupported/Eigen/FFT>

// Classes providing the nodes and weights of Clenshaw-Curtis-Fejer quadrature
// rules on the reference interval [-1,1]
/* SAM_LISTING_BEGIN_1 */
class CCFQuadRule {
public:
  explicit CCFQuadRule(unsigned int n);
  ~CCFQuadRule() = default;

  const Eigen::VectorXd &nodes() const { return nodes_; }
  const Eigen::VectorXd &weights() const { return weights_; }

private:
  Eigen::VectorXd nodes_;
  Eigen::VectorXd weights_;
};
/* SAM_LISTING_END_1 */

// Class with efficient implementation of constructor
/* SAM_LISTING_BEGIN_2 */
class CCFQuadRule_Fast {
public:
  explicit CCFQuadRule_Fast(unsigned int n);
  ~CCFQuadRule_Fast() = default;

  const Eigen::VectorXd &nodes() const { return nodes_; }
  const Eigen::VectorXd &weights() const { return weights_; }

private:
  Eigen::VectorXd nodes_;
  Eigen::VectorXd weights_;
};
/* SAM_LISTING_END_2 */

// Computation of Clenshaw-Curtis Fejer weights and nodes
// in a straightforward way by using Formula (2.3) in WAV06
/* SAM_LISTING_BEGIN_4 */
CCFQuadRule::CCFQuadRule(unsigned int n) : weights_(n + 1) {
  assert(n > 0);
  nodes_ = (Eigen::ArrayXd::LinSpaced(n + 1, 0, n) * M_PI / n).cos().matrix();
  int m = n / 2; // Integer division!
  const Eigen::ArrayXd idx = Eigen::ArrayXd::LinSpaced(m, 1.0, m);
  const Eigen::ArrayXd cos_arg = 2.0 * M_PI * idx / n;
  Eigen::ArrayXd fac = 2.0 / (4 * idx.pow(2) - 1.0);
  if (n % 2 == 0)
    fac[m - 1] /= 2.0;
  weights_[0] = (1.0 - fac.sum()) / n;
  for (int j = 1; j < n; ++j) {
    weights_[j] = (1.0 - (fac * (j * cos_arg).cos()).sum()) * 2.0 / n;
  }
  weights_[n] = weights_[0];
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_3 */
CCFQuadRule_Fast::CCFQuadRule_Fast(unsigned int n) : weights_(n + 1) {
  // TO DO: Provide an implementation with asymptotic complexity of O(n log(n))
  // START
  assert(n > 0);
  // The nodes $c_j = \cos\left(\frac{j\pi}{n}\right)$, $j=0,\ldots,n$
  nodes_ = (Eigen::ArrayXd::LinSpaced(n + 1, 0, n) * M_PI / n).cos().matrix();
  // Length of sum: $m= \frac{n}{2}$ for even $n$, $m = \frac{n-1}{2}$ for odd
  // $n$
  int m = n / 2;
  // Initialize vector $\Vz$ to be multiplied with Fourier matrix
  Eigen::VectorXcd z{Eigen::VectorXcd::Zero(n)};
  z[0] = 1.0;
  for (int l = 1; l <= m; ++l)
    z[l] = -2.0 / (4.0 * l * l - 1.0);
  if (n % 2 == 0)
    z[m] /= 2.0;
  // FFT-based multiplication with Fourier matrix
  Eigen::FFT<double> fft;
  // Eigen's built-in FFT crashes for vectors of length = 1
  // This is a work around
  if (n == 1)
    weights_[0] = 2.0;
  else
    weights_.head(n) = 2.0 * (fft.fwd(z)).real() / n;
  weights_[0] /= 2.0;
  weights_[n] = weights_[0];
  // END
}
/* SAM_LISTING_END_3 */

// Dermination of order of CCF quadrature rule
/* SAM_LISTING_BEGIN_5 */
unsigned int determineOrderCCF(unsigned int n) {
  const CCFQuadRule ccfqr(n);
  int d = 0;       // Degree of monomial, whose integral value is being tested
  double quaderr;  // Quadrature error
  double I_ex;     // Exact integral value
  double I_qr;     // Value computed with CCF quadrature rule
  // Predefined fixed tolerances
  const double abstol = std::numeric_limits<double>::epsilon() * n * 2;
  const double reltol = std::numeric_limits<double>::epsilon() * n * 20;
  // Main quadrature loop
  do {
    // Compute powers of nodes: $\cob{c_j^d}$
    const Eigen::ArrayXd pdvals = ccfqr.nodes().array().pow(d);
    // Evaluate quadrature formula: Implicit summation via dot product
    I_qr = ccfqr.weights().dot(pdvals.matrix());
    // Exact integral value according to \prbeqref{eq:Iex}
    I_ex = (d % 2 == 0) ? (2.0 / (d + 1)) : 0.0;
    quaderr = std::abs(I_qr - I_ex);
    d = d + 1;
  } while ((quaderr <= abstol) || (quaderr <= reltol * I_ex)); // \Label[line]{ccftest}
  return d - 1;
}
/* SAM_LISTING_END_5 */
