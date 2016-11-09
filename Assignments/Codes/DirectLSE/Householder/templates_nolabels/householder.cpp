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
    // TODO
}

template <typename Scalar>
void applyHouseholder(VectorXd &x,
                      const MatrixBase<Scalar> &V) {
    int n = V.rows();

    // TODO
}

int main() {

}
