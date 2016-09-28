//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

int main() {
  auto h = ArrayXd::LinSpaced(21, -20, 0).unaryExpr([] (const double i) {
    return std::pow(10, i);
  });
  auto x = ArrayXd::Constant(h.size(), 1.2);

  // Derivative
  ArrayXd g1 = (sin(x+h) - sin(x)) / h; // naive
  ArrayXd g2 = 2 * cos(x+0.5*h) * sin(0.5 * h) / h;
  ArrayXd ex = cos(x);

  // Plot
  mgl::Figure fig;
  fig.setlog(true, true);
  fig.legend();
  fig.title("Error of approximation of f'(x_0)");
  fig.xlabel("h");
  fig.ylabel("| f'(x_0) - g_i(x_0, h) |");
  fig.plot(h.matrix(), (g1-ex).abs().matrix()).label("g_1");
  fig.plot(h.matrix().tail(16), (g2-ex).abs().matrix().tail(16)).label("g_2");
  fig.plot(h.matrix(), h.matrix(), "h;").label("O(h)");
  fig.save("error_cancellation.eps");
}
