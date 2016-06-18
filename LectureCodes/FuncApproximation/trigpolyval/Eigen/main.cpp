# include "./trigpolyval.hpp"

int main() {
  const unsigned N = 5;
  VectorXd t = VectorXd::LinSpaced(N, 0, 2*M_PI),
           y = t,
           x = VectorXd::LinSpaced(10, 0, 2*M_PI),
           q;
  trigpolyval(t, y, x, q);
  std::cout << q.real() << "\n";
  return 0;
}
