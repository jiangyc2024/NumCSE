#include "meshgrid.hpp" // provided by NumCSE/Utils
#include <Eigen/Dense>

namespace psf {


using Eigen::MatrixXd; 
using Eigen::VectorXd;

inline
/* SAM_LISTING_BEGIN_0 */
void psf(const Eigen::Index L, MatrixXd& S) {
  const VectorXd x = VectorXd::LinSpaced(2*L+1, -static_cast<double>(L), static_cast<double>(L));
  const MatrixXd X = x.replicate(1,x.size());
  const MatrixXd Y = (x.transpose()).replicate(x.size(),1);
  const MatrixXd E = MatrixXd::Ones(2*L+1, 2*L+1);
  S = E.cwiseQuotient(E + X.cwiseProduct(X) + Y.cwiseProduct(Y));
  S /= S.sum();
}
/* SAM_LISTING_END_0 */


} //namespace psf