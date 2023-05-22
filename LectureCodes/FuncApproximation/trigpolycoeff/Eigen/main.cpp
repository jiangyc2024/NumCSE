# include "./trigpolycoeff.hpp"

int main() {
  const Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(11, 0, 0.5);
  const Eigen::VectorXd y = t.cwiseProduct(t);
  Eigen::VectorXd a;
  Eigen::VectorXd b;
  std::tie(a,b) = trigpolycoeff::trigpolycoeff(t, y);
  std::cout << "Alphas: " << a.transpose() << "\n"
            << "Betas: " << b.transpose() << "\n";

  return 0;
}
