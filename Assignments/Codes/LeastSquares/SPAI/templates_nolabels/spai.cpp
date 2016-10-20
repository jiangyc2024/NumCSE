//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <cmath>
#include <ctime>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <unsupported/Eigen/KroneckerProduct>

#include "timer.h"

using namespace Eigen;

using index_t = int;

/* @brief Compute $B = \argmin_{X \in P(A)} |I-AX|_F$
 * @param[in] A An $n \times n$ matrix
 * @param[out] B The $n \times n$ matrix $= \argmin_{X \in P(A)} |I-AX|_F$
 */
SparseMatrix<double> spai(SparseMatrix<double> & A) {
    // Size check
    assert(A.rows() == A.cols() &&
    "Matrix must be square!");
    unsigned int N = A.rows();

    // Needed to make sure ___Ptr functions return
    // arrays specified in CRS format
    A.makeCompressed();

    // Obtain pointers to data of A
    double* valPtr = A.valuePtr();
    index_t* innPtr = A.innerIndexPtr();
    index_t* outPtr = A.outerIndexPtr();

    // Create vector for triplets of B and reserve enough space
    std::vector<Triplet<double>> triplets;
    triplets.reserve(A.nonZeros());

    // TODO: build $B$ in triplet format, exploiting the sparse format of $A$

    // Build and return SPAI preconditioner
    SparseMatrix<double> B = SparseMatrix<double>(N,N);
    B.setFromTriplets(triplets.begin(), triplets.end());
    B.makeCompressed();
    return B;
}

// Run conditionally only one test
const bool small_test = false;
const bool big_test = true;

int main(int argc, char **argv) {
    unsigned int n = 60;
    if(argc > 1) {
        n = std::stoi(argv[1]);
    }
    srand(time(NULL));

    // First test: compute SPAI with very small matrix
    // Check correctness
    if(small_test)
    {
        SparseMatrix<double> M(5,5);
        M.coeffRef(3,4) = 1;
        M.coeffRef(4,4) = 2;
        M.coeffRef(1,4) = 3;
        M.coeffRef(3,3) = 4;
        M.coeffRef(3,2) = 4;
        M.coeffRef(2,3) = 4;
        M.coeffRef(2,2) = 5;
        M.coeffRef(3,1) = 6;
        M.coeffRef(0,0) = 9;

        SparseMatrix<double> N = spai(M);
        SparseMatrix<double>I(n,n);
        I.setIdentity();

        std::cout << "M = " << std::endl
                  << M
                  << std::endl;
        std::cout << "N = " << std::endl
                  << N
                  << std::endl;

        std::cout << "Error:"
                  << (I - M*N).norm()
                  << std::endl;

    }
    // Big test: test with large, sparse matrix
    if(big_test)
    {
        SparseMatrix<double> M(n*n,n*n);

        SparseMatrix<double>I(n,n);
        I.setIdentity();
        MatrixXd R = MatrixXd::Random(n,n);
        M = kroneckerProduct(R, I);

        Timer tm;
        tm.start();
        SparseMatrix<double> N = spai(M);
        tm.stop();

        SparseMatrix<double>Ibig(n*n,n*n);
        Ibig.setIdentity();

        std::cout << "Error (n = " << n*n << "):"
                  << (Ibig - M*N).norm()
                  << std::endl
                  << "Elapsed:"
                  << tm.duration() << " s"
                  << std::endl;
    }
}
