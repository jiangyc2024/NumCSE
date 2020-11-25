///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting gradit.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>

using namespace Eigen;

VectorXd gradit(MatrixXd A, VectorXd b, VectorXd x, double rtol, double atol, unsigned int maxit);