# include "trigipequid.hpp"

int main() {
  const Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(5, 0, 1);
  Eigen::VectorXcd a;
  Eigen::VectorXcd b;
  trigipequid::trigipequid(y, a, b);
  std::cout << "a = " << a.transpose() << "\n"
            << "b = " << b.transpose() << "\n";
  Eigen::VectorXd ar;
  Eigen::VectorXd br;
  std::tie(ar,br) = trigipequid::trigipequid(y);
  std::cout << "Re(a) = " << ar.transpose() << "\n"
	    << "Re(b) = " << br.transpose() << "\n";
  return 0;
}
