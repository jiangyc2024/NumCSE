///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "corzero.hpp"
#include <Eigen/Dense>
#include <iostream>

int main()
{
  Eigen::VectorXd rates;
  Eigen::VectorXd err;

  corzero::fpit(0.4, rates, err);

  const Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(err.size(), 0, static_cast<double>(err.size()-1));

  std::cout << rates << std::endl;
  std::cout << err << std::endl;
  std::cout << x << std::endl;

  return 0;
}
