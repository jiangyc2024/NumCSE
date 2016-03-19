# include <Eigen/Dense>

using Eigen::VectorXd;
// t, y: Interpolation points
// x: Evaluation point
double ANipoleval(const VectorXd& t, VectorXd y, const double& x) {
  for (int i = 0; i < y.size(); ++i) {
    for (int k = i - 1; k >= 0; --k) {
      y(k) = y(k + 1) + (y(k + 1) - y(k))*(x - t(i))/(t(i) - t(k));
    }
  }
  return y(0);
}
