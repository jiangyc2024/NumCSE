///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting pcgbase.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

std::tuple<VectorXd, std::vector<double>, std::vector<VectorXd>> pcgbase(
    std::function<VectorXd(VectorXd)> evalA,
    VectorXd b, double tol, unsigned int maxit,
    std::function<VectorXd(VectorXd)> invB, VectorXd x);