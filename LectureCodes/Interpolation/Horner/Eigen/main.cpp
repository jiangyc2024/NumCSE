# include <iostream>
# include <Eigen/Dense>
# include "./horner.hpp"
# include <chrono>
# include <iomanip>

# include <figure/figure.hpp>

using std::chrono::high_resolution_clock;
# define duration_cast std::chrono::duration_cast<std::chrono::nanoseconds>

void simple_eval(const Eigen::VectorXd& c, const Eigen::VectorXd& x, Eigen::VectorXd& y) {
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
  y = 0*ones;
  for (unsigned i = 0; i < c.size(); ++i) {
    y += c(i)*(x.array().pow(c.size() - 1 - i)).matrix();
  }
}

int main () {

  std::cout << "=========================================================\n"
            << "Testing runtime and values of simple eval and horner eval\n"
            << "=========================================================\n\n";

  std::cout << "**** Number of eval points fixed (50) & coeffs increasing ****\n\n";
  std::cout << "# coeffs     time simple          time horner         difference of results\n";
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(50, 0, 1);


  for (unsigned n = 2; n <= 1000; n *= 2) {
    Eigen::VectorXd c = Eigen::VectorXd::Random(n);
    Eigen::VectorXd ys, yh;
    auto t = high_resolution_clock::now();
    simple_eval(c, x, ys);
    const double ts = duration_cast(high_resolution_clock::now() - t).count()/1e9;
    t = high_resolution_clock::now();
    horner(c, x, yh);
    const double th = duration_cast(high_resolution_clock::now() - t).count()/1e9;

    std::cout << std::setw(4) << n << std::setw(20) << ts << std::setw(20) << th << std::setw(25) << (yh - ys).cwiseAbs().sum() << "\n";
  }

  std::cout << "\n\n**** Number of eval points increasing & coeffs fixed (20) ****\n\n";
  std::cout << "# evals      time simple          time horner         difference of results\n";

  Eigen::VectorXd c = Eigen::VectorXd::Random(20);
  for (unsigned n = 10; n <= 5000; n *= 2) {
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, 0, 1);
    Eigen::VectorXd ys, yh;
    auto t = high_resolution_clock::now();
    simple_eval(c, x, ys);
    const double ts = duration_cast(high_resolution_clock::now() - t).count()/1e9;
    t = high_resolution_clock::now();
    horner(c, x, yh);
    const double th = duration_cast(high_resolution_clock::now() - t).count()/1e9;

    std::cout << std::setw(4) << n << std::setw(20) << ts << std::setw(20) << th << std::setw(25) << (yh - ys).cwiseAbs().sum() << "\n";
  }
  return 0;
}
