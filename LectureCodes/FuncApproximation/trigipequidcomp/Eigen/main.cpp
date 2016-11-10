# include "./../../trigipequid/Eigen/trigipequid.hpp"
# include "./trigipequidcomp.hpp"

int main() {
  const unsigned N = 10;
  VectorXd y = VectorXd::LinSpaced(5, 0, 1);
  VectorXcd a, b, q;
  std::tie(a,b) = trigipequid(y);
  q = trigipequidcomp(a, b, N);

  std::cout << q << "\n";

  return 0;
}
