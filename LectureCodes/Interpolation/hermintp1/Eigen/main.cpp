# include <cmath>
# include <Eigen/Dense>
# include "./hermintp1.hpp"

int main() {
  auto f = [](double x) { return std::sin(5*x)*std::exp(x); };
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(10, -1, 1);

  hermintp1(f, t);

  return 0;
}
