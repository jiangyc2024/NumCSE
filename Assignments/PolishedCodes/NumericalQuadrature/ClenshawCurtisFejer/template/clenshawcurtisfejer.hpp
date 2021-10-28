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
/* SAM_LISTING_BEGIN_1 */
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
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_3 */
CCFQuadRule_Fast::CCFQuadRule_Fast(unsigned int n) : weights_(n + 1) {
  // TO DO: Provide an implementation with asymptotic complexity of O(n log(n))
  // START
  
  // END
}
/* SAM_LISTING_END_3 */
