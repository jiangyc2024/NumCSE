# include <Eigen/Dense>
# include "meshgrid.hpp" // provided by NumCSE/Utils
using Eigen::MatrixXd; using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
void psf(const long L, MatrixXd& S) {
  const VectorXd x = VectorXd::LinSpaced(2*L+1, -L, L);
  const MatrixXd X = x.replicate(1,x.size());
  const MatrixXd Y = (x.transpose()).replicate(x.size(),1);
  MatrixXd E = MatrixXd::Ones(2*L+1, 2*L+1);
  S = E.cwiseQuotient(E + X.cwiseProduct(X) + Y.cwiseProduct(Y));
  S /= S.sum();
}
/* SAM_LISTING_END_0 */
