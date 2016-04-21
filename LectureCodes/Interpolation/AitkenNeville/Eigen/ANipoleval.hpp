# include <Eigen/Dense>

using Eigen::VectorXd;
// IN:  t, y: Interpolation points
//      x: Evaluation point
// OUT: value of interpolant in x
double ANipoleval(const VectorXd& t, VectorXd y, const double& x) {
  for (int i = 0; i < y.size(); ++i) {
    for (int k = i - 1; k >= 0; --k) {
      y(k) = y(k + 1) + (y(k + 1) - y(k))*(x - t(i))/(t(i) - t(k));
    }
  }
  return y(0);
}
