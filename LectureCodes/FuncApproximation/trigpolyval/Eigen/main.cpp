# include "./trigpolyval.hpp"

int main() {
  const unsigned N = 5;
  const Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(N, 0, 2*M_PI);
  const Eigen::VectorXd & y = t;
  const Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(10, 0, 2*M_PI);
  Eigen::VectorXd q;
  Eigen::VectorXd q1;
  trigpolyval::trigpolyval(t, y, x, q);
  q1 = trigpolyval::trigpolyval(t, y, x);
  std::cout << " q = " << q.real().transpose() << "\n";
  std::cout << " q1 = " << q1.real().transpose() << "\n";
  return 0;
}
