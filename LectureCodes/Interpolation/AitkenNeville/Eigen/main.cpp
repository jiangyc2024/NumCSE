# include "./ANipoleval.hpp"
# include <iostream>

int main () {
  // approximate d/dx sin(x) with aitken neville and extrapolation to 0
  const double x = M_PI_2;
  auto df = [x](double h){ return (std::sin(x + h) - std::sin(x))/h; };

  const unsigned N = 50;
  // exact value
  const double exact = std::cos(x);
  std::cout << "Exact value: " << exact << "\n";

  for (unsigned n = 2; n < N; n += 2) {
    // interpolate for large h
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n, 0.01, 10);
    Eigen::VectorXd y(t.size());
    for (unsigned i = 0; i < t.size(); ++i) {
      y(i) = df(t(i));
    }

    std::cout << "Approximation at n = " << n << ": " << ANipoleval::ANipoleval(t, y, 0) << "\n";
  }

  return 0;
}
