# include "./trigpolycoeff.hpp"

int main() {
  VectorXd t = VectorXd::LinSpaced(11, 0, 0.5),
           y = t.cwiseProduct(t),
           a, 
           b;

  trigpolycoeff(t, y, a, b);
  std::cout << "Alphas: " << a.transpose() << "\n"
            << "Betas: " << b.transpose() << "\n";

  return 0;
}
