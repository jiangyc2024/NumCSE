# include "./trigipequid.hpp"

int main() {
  VectorXd y = VectorXd::LinSpaced(5, 0, 1);
  VectorXcd a, b;
  trigipequid(y, a, b);
  std::cout << "a = " << a.transpose() << "\n"
            << "b = " << b.transpose() << "\n";
  VectorXd ar,br;
  std::tie(ar,br) = trigipequid(y);
  std::cout << "Re(a) = " << ar.transpose() << "\n"
	    << "Re(b) = " << br.transpose() << "\n";
  return 0;
}
