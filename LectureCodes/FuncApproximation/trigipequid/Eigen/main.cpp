# include "./trigipequid.hpp"

int main() {

  {
	  // non complex
	  VectorXd y = VectorXd::LinSpaced(5, 0, 1);
	  VectorXcd a, b;
	  std::tie(a,b) = trigipequid(y);
	  std::cout << "a = " << a.transpose() << "\n"
				<< "b = " << b.transpose() << "\n";
  }

  {
	  // complex type
	  VectorXcd y = VectorXd::LinSpaced(5, 0, 1).cast<std::complex<double>>();
	  VectorXcd a, b;
	  std::tie(a,b) = trigipequid(y);
	  std::cout << "a = " << a.transpose() << "\n"
				<< "b = " << b.transpose() << "\n";
  }

  return 0;
}
