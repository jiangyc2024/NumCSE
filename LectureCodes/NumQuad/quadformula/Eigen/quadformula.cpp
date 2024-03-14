///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2020 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "quadformula.hpp"
#include <iostream>

int main(int /*argc*/, char **/*argv*/) {
  const Eigen::VectorXd c { (Eigen::VectorXd(4) << 0.0,0.33,0.66,1.0).finished() };
  const Eigen::VectorXd w { 0.25 * Eigen::VectorXd::Ones(4) };
  std::cout << "integral = " << quadformula([] (double t) { return t*t; }, c, w) << std::endl;
  return 0;
}
