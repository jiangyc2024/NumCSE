# include "./../../trigipequid/Eigen/trigipequid.hpp"
# include "./trigipequidcomp.hpp"

int main() {
  const unsigned N = 10;
  const Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(5, 0, 1);
  Eigen::VectorXcd a;
  Eigen::VectorXcd b;
  Eigen::VectorXcd q;
  trigipequid::trigipequid(y, a, b);
  trigipequidcomp::trigipequidcomp(a, b, N, q);

  std::cout << q << "\n";

  return 0;
}
