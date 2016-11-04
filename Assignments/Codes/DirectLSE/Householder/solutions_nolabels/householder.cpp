//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

void applyHouseholder(VectorXd &x, const VectorXd &v) {
    double d = v.transpose()*x;
    x -= 2*v*d;
}

template <typename Scalar>
void applyHouseholder(VectorXd &x,
                      const MatrixBase<Scalar> &V) {
    int n = V.rows();

    for(unsigned int i = 0; i < n; ++i) {
        VectorXd v = V.col(i);
        v(i) = std::sqrt(1. - V.col(i).normSquared());

        x += 2.*v.dot(x) / (1. + 2.*v.dot(v)) * v;
    }
}

int main() {

}
