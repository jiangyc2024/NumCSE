#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;


void applyHouseholder(VectorXd &x, const VectorXd &v) {
    double d = v.transposed()*x;
    x -= 2*v*d;
}

template <typename Scalar>
void applyHouseholder(VectorXd &x,
                      const MatrixBase<Scalar> &V) {

}



int main() {

}
