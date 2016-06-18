# include "./trigpolyvalequid.hpp"

int main() {
  VectorXd y = VectorXd::LinSpaced(11, 0, 1).array().pow(2).matrix();
  VectorXcd q;
  trigpolyvalequid(y, q);
  std::cout << q << "\n";

  return 0;
}
