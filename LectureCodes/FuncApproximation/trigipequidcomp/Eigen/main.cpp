# include "./../../trigipequid/Eigen/trigipequid.hpp"
# include "./trigipequidcomp.hpp"

int main() {
  const unsigned N = 10;
  VectorXd y = VectorXd::LinSpaced(5, 0, 1);
  VectorXcd a, b, q;
  trigipequid(y, a, b);
  trigipequidcomp(a, b, N, q);

  std::cout << q << "\n";

  return 0;
}
