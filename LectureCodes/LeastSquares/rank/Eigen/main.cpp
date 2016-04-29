#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main()
{
    MatrixXd A(3, 2);
    double eps = NumTraits<double>::epsilon();
    A << 1, 1,
         sqrt(eps), 0,
         0, sqrt(eps);

     cout << A.fullPivLu().rank() << endl;
     cout << (A.transpose() * A).fullPivLu().rank() << endl;

    return 0;
}
