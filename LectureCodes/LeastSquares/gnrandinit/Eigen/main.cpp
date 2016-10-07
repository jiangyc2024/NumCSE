#include "gnrandinit.hpp"
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main()
{
    Eigen::Vector3d x(1., 2., 1.);

    cout << gnrandinit(x) << endl;

    return 0;
}
