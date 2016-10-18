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

SparseMatrix<double> spai(SparseMatrix<double> & A) {
    assert(A.rows() = A.cols() &&
           "Matrix must be square!");

    unsigned int n = A.rows();

    A.makeCompressed();

    double* valPtr = A.valuePtr();
    index_t* innPtr = A.innerIndexPtr();
    index_t* outPtr = A.outerIndexPtr();

    std::vector<Triplet<double>> triplets;
    triplets.reserve(A.nonZeros());

    // Without this, the code is MUCH slower w.r.t. MATLAB
    #pragma omp parallel for
    for(unsigned int i = 0; i < n; ++i) {

        index_t nnz_i = outPtr[i+1] - outPtr[i];

        if(nnz_i == 0) continue;

        MatrixXd C = MatrixXd::Zero(n, nnz_i);
        for(unsigned int k = outPtr[i]; k < outPtr[i+1]; ++k) {
            index_t row_k = innPtr[k];
            index_t nnz_k = outPtr[row_k+1] - outPtr[row_k];
            for(unsigned int l = 0; l < nnz_k; ++l) {
                unsigned int innIdx = outPtr[row_k] + l;
                C(innPtr[innIdx], k - outPtr[i]) = valPtr[innIdx];
            }
        }
        VectorXd b = (C.transpose() * C).partialPivLu().solve(C.row(i).transpose());
        for(unsigned int k = 0; k < b.size(); ++k) {
            #pragma omp critical(triplets)
            {
				triplets.push_back(Triplet<double>(innPtr[outPtr[i] + k], i, b(k)));
            }
        }

    }

    SparseMatrix<double> B = SparseMatrix<double>(n,n);
    B.setFromTriplets(triplets.begin(), triplets.end());
    B.makeCompressed();

//    std::cout << B;

    return B;
}

// Run conditionally
const bool small_test = false;
const bool big_test = true;

int main(int argc, char **argv) {
    unsigned int n = 60;
    if(argc > 1) {
        n = std::stoi(argv[1]);
    }
    srand(time(NULL));

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
    if(big_test)
    {
        SparseMatrix<double> M(n*n,n*n);

        SparseMatrix<double>I(n,n);
        I.setIdentity();
        MatrixXd R = MatrixXd::Random(n,n);
//        M = kroneckerProduct(I, R);
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
