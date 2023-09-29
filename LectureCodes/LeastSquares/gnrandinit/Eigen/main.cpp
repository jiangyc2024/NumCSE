///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Xiaolin Guo, Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "gnrandinit.hpp"
#include <Eigen/Dense>
#include <iostream>

using std::cout;
using std::endl;

int main()
{
    const Eigen::Vector3d x(1., 2., 1.);

    cout << gnrandinit::gnrandinit(x) << endl;

    return 0;
}
