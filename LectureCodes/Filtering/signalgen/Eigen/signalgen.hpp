# include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::ArrayXd;

/* SAM_LISTING_BEGIN_0 */
VectorXd signalgen() {
  const int N = 64;
  const ArrayXd t = ArrayXd::LinSpaced(N, 0, N);
  const VectorXd x = ((2*M_PI/N*t).sin() + (14*M_PI/N*t).sin()).matrix();
  return x + VectorXd::Random(N);
}
/* SAM_LISTING_END_0 */
