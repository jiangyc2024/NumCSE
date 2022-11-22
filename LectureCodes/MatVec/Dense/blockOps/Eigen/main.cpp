///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>

#include "./blockOps.hpp"

int main () {
  Eigen::MatrixXd M = (10*Eigen::MatrixXd::Random(4,4)).cast<int>().cast<double>();
  blockops::blockAccess(M);
  return 0;
}
