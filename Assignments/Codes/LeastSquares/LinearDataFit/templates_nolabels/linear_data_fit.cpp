//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <cmath>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

MatrixXd make_A(const VectorXd &b) {
	size_t n = b.size();
	MatrixXd A(n, 4);
	// TODO: construct fitting matrix A
	return A;
}


VectorXd data_fit_normal(const MatrixXd &A, const VectorXd &b) {
	// TODO: solve least-squares problem using multiplication with transpose
}

VectorXd data_fit_qr(const MatrixXd &A, const VectorXd &b) {
	// TODO: solve least-squares problem using QR-Decomposition
}


int main() {
	// TODO: create a lin-log plot using the mgl::Figure class

	// TODO: analyze errors
}
