#include <Eigen/Core>
#include <iostream>

#include "linearinterpolant.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
  // Test the class with the basis with nodes (2,-1,4,1)
  // and interpolant with values (3,-1,4,2).
  Eigen::VectorXd t(4);
  t << 2., -1., 4., 1.;
  Eigen::VectorXd y(4);
  y << 3., -1., 4., 2.;
  LinearInterpolant I(t, y);

  std::cout << "Testing with (t_i, y_i) = (2, 3), (-1, -1), (4, 4), (1, 2).\n";

  int n = 100;
  Eigen::VectorXd x_values = Eigen::VectorXd::LinSpaced(n, -1, 4);
  Eigen::VectorXd y_interpolated =
      x_values.unaryExpr([&I](double x) { return I(x); });

  plt::figure();
  plt::title("Piecewise linear interpolation polynomial.");
  plt::plot(t, y, "r*", {{"label", "data points"}});
  plt::plot(x_values, y_interpolated, "b", {{"label", "I(x)"}});
  plt::xlabel("x");
  plt::ylabel("y");
  plt::legend();
  plt::savefig("cx_out/linearinterpolant.png");

  std::cout << "A plot has been generated.\n";
}
