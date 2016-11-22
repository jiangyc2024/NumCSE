# include "./trigpolyval.hpp"

int main() {
  const unsigned N = 5;
  VectorXd t = VectorXd::LinSpaced(N, 0, 2*M_PI),
    y = t,
    x = VectorXd::LinSpaced(10, 0, 2*M_PI),
    q,q1;
  trigpolyval(t, y, x, q);
  q1 = trigpolyval(t, y, x);
  std::cout << " q = " << q.real().transpose() << "\n";
  std::cout << " q1 = " << q1.real().transpose() << "\n";
  return 0;
}
