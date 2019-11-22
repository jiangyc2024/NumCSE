#include <Eigen/Dense>
using Eigen::VectorXd;

VectorXd gnrandinit(const VectorXd &x) {
  std::srand((unsigned int)time(0));
  auto t = VectorXd::LinSpaced((7.0 - 1.0) / 0.3 - 1, 1.0, 7.0);
  auto y = x(0) + x(1) * ((-x(2) * t).array().exp());
  return y + 0.1 * (VectorXd::Random(y.size()).array() - 0.5);
}
