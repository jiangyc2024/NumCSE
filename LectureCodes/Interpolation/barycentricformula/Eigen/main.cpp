///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): J. Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <figure/figure.hpp>
#include <iostream>

#include "./intpolyval.hpp"

int main() {
  VectorXd t = VectorXd::LinSpaced(10, -5, 5),
           y = (1. / (1 + t.array() * t.array())).matrix();

  VectorXd x = VectorXd::LinSpaced(100, -5, 5), p;

  p = intpolyval(t, y, x);

  mgl::Figure fig;
  fig.plot(t, y, " *n").label("Data");
  fig.plot(x, p, "g").label("Interpolant");
  fig.fplot("1/(1 + x^2)", "r").label("Function");
  fig.legend();
  fig.save("interpolation");

  return 0;
}
