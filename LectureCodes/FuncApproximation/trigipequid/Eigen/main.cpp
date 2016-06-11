# include "./trigipequid.hpp"

int main() {
  VectorXd y = VectorXd::LinSpaced(5, 0, 1);
  VectorXcd a, b;
  trigipequid(y, a, b);
  std::cout << "a = " << a.transpose() << "\n"
            << "b = " << b.transpose() << "\n";

  return 0;
}
